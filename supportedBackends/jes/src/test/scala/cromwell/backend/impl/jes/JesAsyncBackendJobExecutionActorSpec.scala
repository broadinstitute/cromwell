package cromwell.backend.impl.jes

import java.util.UUID

import _root_.io.grpc.Status
import _root_.wdl._
import wdl.transforms.draft2.wdlom2wom._
import wom.transforms.WomExecutableMaker.ops._
import wom.transforms.WomWorkflowDefinitionMaker.ops._
import akka.actor.{ActorRef, Props}
import akka.testkit.{ImplicitSender, TestActorRef, TestDuration, TestProbe}
import com.google.api.gax.retrying.RetrySettings
import com.google.cloud.NoCredentials
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, JobFailedNonRetryableResponse, JobFailedRetryableResponse}
import cromwell.backend._
import cromwell.backend.async.AsyncBackendJobExecutionActor.{Execute, ExecutionMode}
import cromwell.backend.async.{AbortedExecutionHandle, ExecutionHandle, FailedNonRetryableExecutionHandle, FailedRetryableExecutionHandle}
import cromwell.backend.impl.jes.JesAsyncBackendJobExecutionActor.JesPendingExecutionHandle
import cromwell.backend.impl.jes.RunStatus.UnsuccessfulRunStatus
import cromwell.backend.impl.jes.io.{DiskType, JesWorkingDisk}
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager.DoPoll
import cromwell.backend.standard.{DefaultStandardAsyncExecutionActorParams, StandardAsyncExecutionActorParams, StandardAsyncJob, StandardExpressionFunctionsParams}
import cromwell.cloudsupport.gcp.gcs.GcsStorage
import cromwell.core.Tags.PostWomTest
import cromwell.core._
import cromwell.core.callcaching.NoDocker
import cromwell.core.labels.Labels
import cromwell.core.logging.JobLogger
import cromwell.core.path.{DefaultPathBuilder, PathBuilder}
import cromwell.filesystems.gcs.{GcsPath, GcsPathBuilder}
import cromwell.services.keyvalue.InMemoryKvServiceActor
import cromwell.services.keyvalue.KeyValueServiceActor.{KvJobKey, KvPair, ScopedKey}
import cromwell.util.JsonFormatting.WomValueJsonFormatter._
import cromwell.util.SampleWdl
import org.scalatest._
import org.scalatest.prop.Tables.Table
import org.slf4j.Logger
import org.specs2.mock.Mockito
import spray.json._
import wom.WomFileMapper
import wom.graph.CommandCallNode
import wom.types._
import wom.values._

import scala.concurrent.duration._
import scala.concurrent.{Await, ExecutionContext, Future, Promise}
import scala.language.postfixOps
import scala.util.Success

class JesAsyncBackendJobExecutionActorSpec extends TestKitSuite("JesAsyncBackendJobExecutionActorSpec")
  with FlatSpecLike with Matchers with ImplicitSender with Mockito with BackendSpec with BeforeAndAfter with DefaultJsonProtocol {
  val mockPathBuilder: GcsPathBuilder = GcsPathBuilder.fromCredentials(NoCredentials.getInstance(),
    "test-cromwell", RetrySettings.newBuilder().build(), GcsStorage.DefaultCloudStorageConfiguration, WorkflowOptions.empty)
  
  var kvService: ActorRef = system.actorOf(Props(new InMemoryKvServiceActor))

  import JesTestConfig._

  implicit val Timeout: FiniteDuration = 25.seconds.dilated

  val YoSup: String =
    s"""
      |task sup {
      |  String addressee
      |  command {
      |    echo "yo sup $${addressee}!"
      |  }
      |  output {
      |    String salutation = read_string(stdout())
      |  }
      |  runtime {
      |    docker: "ubuntu:latest"
      |    [PREEMPTIBLE]
      |  }
      |}
      |
      |workflow wf_sup {
      |  call sup
      |}
    """.stripMargin

  val Inputs: Map[FullyQualifiedName, WomValue] = Map("wf_sup.sup.addressee" -> WomString("dog"))

  val NoOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue]))

  lazy val TestableCallContext = CallContext(mockPathBuilder.build("gs://root").get, "out", "err")

  lazy val TestableStandardExpressionFunctionsParams = new StandardExpressionFunctionsParams {
    override lazy val pathBuilders: List[PathBuilder] = List(mockPathBuilder)
    override lazy val callContext: CallContext = TestableCallContext
    override val ioActorProxy: ActorRef = simpleIoActor
    override val executionContext = system.dispatcher
  }

  lazy val TestableJesExpressionFunctions: JesExpressionFunctions = {
    new JesExpressionFunctions(TestableStandardExpressionFunctionsParams)
  }

  private def buildInitializationData(jobDescriptor: BackendJobDescriptor, configuration: JesConfiguration) = {
    val workflowPaths = JesWorkflowPaths(jobDescriptor.workflowDescriptor, NoCredentials.getInstance(), NoCredentials.getInstance(), configuration)(system)
    val runtimeAttributesBuilder = JesRuntimeAttributes.runtimeAttributesBuilder(configuration)
    JesBackendInitializationData(workflowPaths, runtimeAttributesBuilder, configuration, null, null)
  }

  class TestableJesJobExecutionActor(params: StandardAsyncExecutionActorParams, functions: JesExpressionFunctions)
    extends JesAsyncBackendJobExecutionActor(params) {

    def this(jobDescriptor: BackendJobDescriptor,
             promise: Promise[BackendJobExecutionResponse],
             jesConfiguration: JesConfiguration,
             functions: JesExpressionFunctions = TestableJesExpressionFunctions,
             jesSingletonActor: ActorRef = emptyActor,
             ioActor: ActorRef = mockIoActor) = {

      this(
        DefaultStandardAsyncExecutionActorParams(
          jobIdKey = JesAsyncBackendJobExecutionActor.JesOperationIdKey,
          serviceRegistryActor = kvService,
          ioActor = ioActor,
          jobDescriptor = jobDescriptor,
          configurationDescriptor = jesConfiguration.configurationDescriptor,
          backendInitializationDataOption = Option(buildInitializationData(jobDescriptor, jesConfiguration)),
          backendSingletonActorOption = Option(jesSingletonActor),
          completionPromise = promise,
          minimumRuntimeSettings = MinimumRuntimeSettings()
        ),
        functions
      )
    }

    override lazy val jobLogger = new JobLogger("TestLogger", workflowId, jobTag, akkaLogger = Option(log)) {
      override def tag: String = s"$name [UUID(${workflowId.shortString})$jobTag]"
      override val slf4jLoggers: Set[Logger] = Set.empty
    }

    override lazy val backendEngineFunctions: JesExpressionFunctions = functions
  }

  private val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)
  private val runtimeAttributesBuilder = JesRuntimeAttributes.runtimeAttributesBuilder(jesConfiguration)
  private val workingDisk = JesWorkingDisk(DiskType.SSD, 200)

  val DockerAndDiskRuntime: String =
    """
      |runtime {
      |  docker: "ubuntu:latest"
      |  disks: "local-disk 200 SSD"
      |}
    """.stripMargin

  private def buildPreemptibleJobDescriptor(preemptible: Int, previousPreemptions: Int, previousUnexpectedRetries: Int): BackendJobDescriptor = {
    val attempt = previousPreemptions + previousUnexpectedRetries + 1
    val wdlNamespace = WdlNamespaceWithWorkflow.load(YoSup.replace("[PREEMPTIBLE]", s"preemptible: $preemptible"),
      Seq.empty[ImportResolver]).get
    val womDefinition = wdlNamespace.workflow.toWomWorkflowDefinition.getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))

    wdlNamespace.toWomExecutable(Option(Inputs.toJson.compactPrint)) match {
      case Right(womExecutable) =>
        val inputs = for {
          combined <- womExecutable.resolvedExecutableInputs
          (port, resolvedInput) = combined
          value <- resolvedInput.select[WomValue]
        } yield port -> value

        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womDefinition,
          inputs,
          NoOptions,
          Labels.empty
        )

        val job = workflowDescriptor.callable.taskCallNodes.head
        val key = BackendJobDescriptorKey(job, None, attempt)
        val runtimeAttributes = makeRuntimeAttributes(job)
        val prefetchedKvEntries = Map(
          JesBackendLifecycleActorFactory.preemptionCountKey -> KvPair(ScopedKey(workflowDescriptor.id, KvJobKey(key), JesBackendLifecycleActorFactory.preemptionCountKey), Some(previousPreemptions.toString)),
          JesBackendLifecycleActorFactory.unexpectedRetryCountKey -> KvPair(ScopedKey(workflowDescriptor.id, KvJobKey(key), JesBackendLifecycleActorFactory.unexpectedRetryCountKey), Some(previousUnexpectedRetries.toString)))
        BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, fqnWdlMapToDeclarationMap(Inputs), NoDocker, prefetchedKvEntries)
      case Left(badtimes) => fail(badtimes.toList.mkString(", "))
    }
  }

  private def executionActor(jobDescriptor: BackendJobDescriptor,
                             configurationDescriptor: BackendConfigurationDescriptor,
                             promise: Promise[BackendJobExecutionResponse],
                             jesSingletonActor: ActorRef,
                             shouldBePreemptible: Boolean): ActorRef = {

    val job = StandardAsyncJob(UUID.randomUUID().toString)
    val run = Run(job, null)
    val handle = new JesPendingExecutionHandle(jobDescriptor, run.job, Option(run), None)

    class ExecuteOrRecoverActor extends TestableJesJobExecutionActor(jobDescriptor, promise, jesConfiguration, jesSingletonActor = jesSingletonActor) {
      override def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
        if(preemptible == shouldBePreemptible) Future.successful(handle)
        else Future.failed(new Exception(s"Test expected preemptible to be $shouldBePreemptible but got $preemptible"))
      }
    }

    system.actorOf(Props(new ExecuteOrRecoverActor), "ExecuteOrRecoverActor-" + UUID.randomUUID)
  }

  private def runAndFail(previousPreemptions: Int, previousUnexpectedRetries: Int, preemptible: Int, errorCode: Status, innerErrorMessage: String, expectPreemptible: Boolean): BackendJobExecutionResponse = {

    val runStatus = UnsuccessfulRunStatus(errorCode, Option(innerErrorMessage), Seq.empty, Option("fakeMachine"), Option("fakeZone"), Option("fakeInstance"), expectPreemptible)
    val statusPoller = TestProbe()

    val promise = Promise[BackendJobExecutionResponse]()
    val jobDescriptor =  buildPreemptibleJobDescriptor(preemptible, previousPreemptions, previousUnexpectedRetries)

    // TODO: Use this to check the new KV entries are there!
    //val kvProbe = TestProbe()

    val backend = executionActor(jobDescriptor, JesBackendConfigurationDescriptor, promise, statusPoller.ref, expectPreemptible)
    backend ! Execute
    statusPoller.expectMsgPF(max = Timeout, hint = "awaiting status poll") {
      case _: DoPoll => backend ! runStatus
    }

    Await.result(promise.future, Timeout)
  }

  def buildPreemptibleTestActorRef(attempt: Int, preemptible: Int): TestActorRef[TestableJesJobExecutionActor] = {
    // For this test we say that all previous attempts were preempted:
    val jobDescriptor = buildPreemptibleJobDescriptor(preemptible, attempt - 1, 0)
    val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(),
      jesConfiguration,
      TestableJesExpressionFunctions,
      emptyActor,
      failIoActor))
    TestActorRef(props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")
  }

  behavior of "JesAsyncBackendJobExecutionActor"

  val timeout = 25 seconds

  { // Set of "handle call failures appropriately with respect to preemption and failure" tests
    val expectations = Table(
      ("previous_preemptions", "previous_unexpectedRetries", "preemptible", "errorCode", "message", "shouldRunAsPreemptible", "shouldRetry"),
      // No preemptible attempts allowed, but standard failures should be retried.
      (0, 0, 0, Status.ABORTED, "13: retryable error", false, true), // This is the new "unexpected failure" mode, which is now retried
      (0, 1, 0, Status.ABORTED, "13: retryable error", false, true),
      (0, 2, 0, Status.ABORTED, "13: retryable error", false, false), // The third unexpected failure is a real failure.
      (0, 0, 0, Status.ABORTED, "14: usually means preempted...?", false, false), // Usually means "preempted', but this wasn't a preemptible VM, so this should just be a failure.
      (0, 0, 0, Status.ABORTED, "15: other error", false, false),
      (0, 0, 0, Status.OUT_OF_RANGE, "13: unexpected error", false, false),
      (0, 0, 0, Status.OUT_OF_RANGE, "14: test error msg", false, false),
      // These commented out tests should be uncommented if/when we stop mapping 13 to 14 in preemption mode
      // 1 preemptible attempt allowed, but not all failures represent preemptions.
//      (0, 0, 1, Status.ABORTED, "13: retryable error", true, true),
//      (0, 1, 1, Status.ABORTED, "13: retryable error", true, true),
//      (0, 2, 1, Status.ABORTED, "13: retryable error", true, false),
      // The following 13 based test should be removed if/when we stop mapping 13 to 14 in preemption mode
      (0, 0, 1, Status.ABORTED, "13: retryable error", true, true),
      (0, 0, 1, Status.ABORTED, "14: preempted", true, true),
      (0, 0, 1, Status.UNKNOWN, "Instance failed to start due to preemption.", true, true),
      (0, 0, 1, Status.ABORTED, "15: other error", true, false),
      (0, 0, 1, Status.OUT_OF_RANGE, "13: retryable error", true, false),
      (0, 0, 1, Status.OUT_OF_RANGE, "14: preempted", true, false),
      (0, 0, 1, Status.OUT_OF_RANGE, "Instance failed to start due to preemption.", true, false),
      // 1 preemptible attempt allowed, but since we're now on the second preemption attempt only 13s should be retryable.
      (1, 0, 1, Status.ABORTED, "13: retryable error", false, true),
      (1, 1, 1, Status.ABORTED, "13: retryable error", false, true),
      (1, 2, 1, Status.ABORTED, "13: retryable error", false, false),
      (1, 0, 1, Status.ABORTED, "14: preempted", false, false),
      (1, 0, 1, Status.UNKNOWN, "Instance failed to start due to preemption.", false, false),
      (1, 0, 1, Status.ABORTED, "15: other error", false, false),
      (1, 0, 1, Status.OUT_OF_RANGE, "13: retryable error", false, false),
      (1, 0, 1, Status.OUT_OF_RANGE, "14: preempted", false, false),
      (1, 0, 1, Status.OUT_OF_RANGE, "Instance failed to start due to preemption.", false, false)
    )

    expectations foreach { case (previousPreemptions, previousUnexpectedRetries, preemptible, errorCode, innerErrorMessage, shouldBePreemptible, shouldRetry) =>
      val descriptor = s"previousPreemptions=$previousPreemptions, previousUnexpectedRetries=$previousUnexpectedRetries preemptible=$preemptible, errorCode=$errorCode, innerErrorMessage=$innerErrorMessage"
      it should s"handle call failures appropriately with respect to preemption and failure ($descriptor)" in {
        runAndFail(previousPreemptions, previousUnexpectedRetries, preemptible, errorCode, innerErrorMessage, shouldBePreemptible) match {
          case response: JobFailedNonRetryableResponse =>
            if(shouldRetry) fail(s"A should-be-retried job ($descriptor) was sent back to the engine with: $response")
          case response: JobFailedRetryableResponse =>
            if(!shouldRetry) fail(s"A shouldn't-be-retried job ($descriptor) was sent back to the engine with $response")
          case huh => fail(s"Unexpected response: $huh")
        }
      }
    }
  }

  it should "not restart 2 of 1 unexpected shutdowns without another preemptible VM" in {
    val actorRef = buildPreemptibleTestActorRef(2, 1)
    val jesBackend = actorRef.underlyingActor
    val runId = StandardAsyncJob(UUID.randomUUID().toString)
    val handle = new JesPendingExecutionHandle(null, runId, None, None)

    val failedStatus = UnsuccessfulRunStatus(Status.ABORTED, Option("14: VM XXX shut down unexpectedly."), Seq.empty, Option("fakeMachine"), Option("fakeZone"), Option("fakeInstance"), true)
    val executionResult = jesBackend.handleExecutionResult(failedStatus, handle)
    val result = Await.result(executionResult, timeout)
    result.isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    val failedHandle = result.asInstanceOf[FailedNonRetryableExecutionHandle]
    failedHandle.returnCode shouldBe None
  }

  it should "restart 1 of 1 unexpected shutdowns without another preemptible VM" in {
    val actorRef = buildPreemptibleTestActorRef(1, 1)
    val jesBackend = actorRef.underlyingActor
    val runId = StandardAsyncJob(UUID.randomUUID().toString)
    val handle = new JesPendingExecutionHandle(null, runId, None, None)

    val failedStatus = UnsuccessfulRunStatus(Status.ABORTED, Option("14: VM XXX shut down unexpectedly."), Seq.empty, Option("fakeMachine"), Option("fakeZone"), Option("fakeInstance"), true)
    val executionResult = jesBackend.handleExecutionResult(failedStatus, handle)
    val result = Await.result(executionResult, timeout)
    result.isInstanceOf[FailedRetryableExecutionHandle] shouldBe true
    val retryableHandle = result.asInstanceOf[FailedRetryableExecutionHandle]
    retryableHandle.returnCode shouldBe None
    retryableHandle.throwable.getMessage should include("will be restarted with a non-preemptible VM")
  }

  it should "restart 1 of 2 unexpected shutdowns with another preemptible VM" in {
    val actorRef = buildPreemptibleTestActorRef(1, 2)
    val jesBackend = actorRef.underlyingActor
    val runId = StandardAsyncJob(UUID.randomUUID().toString)
    val handle = new JesPendingExecutionHandle(null, runId, None, None)

    val failedStatus = UnsuccessfulRunStatus(Status.ABORTED, Option("14: VM XXX shut down unexpectedly."), Seq.empty, Option("fakeMachine"), Option("fakeZone"), Option("fakeInstance"), true)
    val executionResult = jesBackend.handleExecutionResult(failedStatus, handle)
    val result = Await.result(executionResult, timeout)
    result.isInstanceOf[FailedRetryableExecutionHandle] shouldBe true
    val retryableHandle = result.asInstanceOf[FailedRetryableExecutionHandle]
    retryableHandle.returnCode shouldBe None
    retryableHandle.throwable.getMessage should include("will be restarted with another preemptible VM")
  }

  it should "treat a JES message 13 as preemptible if the VM was preemptible" in {
    val actorRef = buildPreemptibleTestActorRef(1, 2)
    val jesBackend = actorRef.underlyingActor
    val runId = StandardAsyncJob(UUID.randomUUID().toString)
    val handle = new JesPendingExecutionHandle(null, runId, None, None)

    val failedStatus = UnsuccessfulRunStatus(Status.ABORTED, Option("13: Retryable Error."), Seq.empty, Option("fakeMachine"), Option("fakeZone"), Option("fakeInstance"), true)
    val executionResult = jesBackend.handleExecutionResult(failedStatus, handle)
    val result = Await.result(executionResult, timeout)
    result.isInstanceOf[FailedRetryableExecutionHandle] shouldBe true
    val retryableHandle = result.asInstanceOf[FailedRetryableExecutionHandle]
    retryableHandle.returnCode shouldBe None
    retryableHandle.throwable.getMessage should include("will be restarted with another preemptible VM")
  }

  it should "handle Failure Status for various errors" in {
    val actorRef = buildPreemptibleTestActorRef(1, 1)
    val jesBackend = actorRef.underlyingActor
    val runId = StandardAsyncJob(UUID.randomUUID().toString)
    val handle = new JesPendingExecutionHandle(null, runId, None, None)

    def checkFailedResult(errorCode: Status, errorMessage: Option[String]): ExecutionHandle = {
      val failed = UnsuccessfulRunStatus(errorCode, errorMessage, Seq.empty, Option("fakeMachine"), Option("fakeZone"), Option("fakeInstance"), true)
      Await.result(jesBackend.handleExecutionResult(failed, handle), timeout)
    }

    checkFailedResult(Status.ABORTED, Option("15: Other type of error."))
      .isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    checkFailedResult(Status.OUT_OF_RANGE, Option("14: Wrong errorCode.")).isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    checkFailedResult(Status.ABORTED, Option("Weird error message.")).isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    checkFailedResult(Status.ABORTED, Option("UnparsableInt: Even weirder error message."))
      .isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    checkFailedResult(Status.ABORTED, None).isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    checkFailedResult(Status.CANCELLED, Option("Operation canceled at")) shouldBe AbortedExecutionHandle

    actorRef.stop()
  }

  it should "map GCS paths and *only* GCS paths to local" taggedAs PostWomTest ignore {
    val stringKey = "abc"
    val stringVal = WomString("abc")
    val localFileKey = "lf"
    val localFileVal = WomSingleFile("/blah/abc")
    val gcsFileKey = "gcsf"
    val gcsFileVal = WomSingleFile("gs://blah/abc")

    val inputs: Map[String, WomValue] = Map(
      stringKey -> stringVal,
      localFileKey -> localFileVal,
      gcsFileKey -> gcsFileVal
    )

    val wdlNamespace = WdlNamespaceWithWorkflow.load(YoSup.replace("[PREEMPTIBLE]", ""),
      Seq.empty[ImportResolver]).get
    val womWorkflow = wdlNamespace.workflow.toWomWorkflowDefinition.getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    wdlNamespace.toWomExecutable(Option(inputs.toJson.compactPrint)) match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap({case (port, v) => v.select[WomValue] map { port -> _ }})

        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty
        )

        val call: CommandCallNode = workflowDescriptor.callable.graph.nodes.collectFirst({ case t: CommandCallNode => t }).get
        val key = BackendJobDescriptorKey(call, None, 1)
        val runtimeAttributes = makeRuntimeAttributes(call)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, fqnWdlMapToDeclarationMap(inputs), NoDocker, Map.empty)

        val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
        val testActorRef = TestActorRef[TestableJesJobExecutionActor](
          props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")


        def gcsPathToLocal(womValue: WomValue): WomValue = {
          WomFileMapper.mapWomFiles(testActorRef.underlyingActor.mapCommandLineWomFile)(womValue).get
        }

        val mappedInputs = jobDescriptor.localInputs mapValues gcsPathToLocal

        mappedInputs(stringKey) match {
          case WomString(v) => assert(v.equalsIgnoreCase(stringVal.value))
          case _ => fail("test setup error")
        }

        mappedInputs(localFileKey) match {
          case wdlFile: WomSingleFile => assert(wdlFile.value.equalsIgnoreCase(localFileVal.value))
          case _ => fail("test setup error")
        }

        mappedInputs(gcsFileKey) match {
          case wdlFile: WomSingleFile => assert(wdlFile.value.equalsIgnoreCase("/cromwell_root/blah/abc"))
          case _ => fail("test setup error")
        }
      case Left(badtimes) => fail(badtimes.toList.mkString(", "))
    }
  }

  private val dockerAndDiskWdlNamespace = WdlNamespaceWithWorkflow.load(SampleWdl.CurrentDirectory.asWorkflowSources(DockerAndDiskRuntime).workflowSource,
    Seq.empty[ImportResolver]).get
  
  it should "generate correct JesFileInputs from a WdlMap" taggedAs PostWomTest ignore {
    val inputs: Map[String, WomValue] = Map(
      "stringToFileMap" -> WomMap(WomMapType(WomStringType, WomSingleFileType), Map(
        WomString("stringTofile1") -> WomSingleFile("gs://path/to/stringTofile1"),
        WomString("stringTofile2") -> WomSingleFile("gs://path/to/stringTofile2")
      )),
      "fileToStringMap" -> WomMap(WomMapType(WomSingleFileType, WomStringType), Map(
        WomSingleFile("gs://path/to/fileToString1") -> WomString("fileToString1"),
        WomSingleFile("gs://path/to/fileToString2") -> WomString("fileToString2")
      )),
      "fileToFileMap" -> WomMap(WomMapType(WomSingleFileType, WomSingleFileType), Map(
        WomSingleFile("gs://path/to/fileToFile1Key") -> WomSingleFile("gs://path/to/fileToFile1Value"),
        WomSingleFile("gs://path/to/fileToFile2Key") -> WomSingleFile("gs://path/to/fileToFile2Value")
      )),
      "stringToString" -> WomMap(WomMapType(WomStringType, WomStringType), Map(
        WomString("stringToString1") -> WomString("path/to/stringToString1"),
        WomString("stringToString2") -> WomString("path/to/stringToString2")
      ))
    )
    
    val workflowInputs = inputs map {
      case (k, v) => s"wf_whereami.whereami$k" -> v
    }

    val womWorkflow = dockerAndDiskWdlNamespace.workflow.toWomWorkflowDefinition.getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    dockerAndDiskWdlNamespace.toWomExecutable(Option(workflowInputs.toJson.compactPrint)) match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap({case (port, v) => v.select[WomValue] map { port -> _ }})
        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty
        )

        val job: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
        val runtimeAttributes = makeRuntimeAttributes(job)
        val key = BackendJobDescriptorKey(job, None, 1)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, fqnWdlMapToDeclarationMap(inputs), NoDocker, Map.empty)

        val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
        val testActorRef = TestActorRef[TestableJesJobExecutionActor](
          props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

        val jesInputs = testActorRef.underlyingActor.generateJesInputs(jobDescriptor)
        jesInputs should have size 8
        jesInputs should contain(JesFileInput(
          "stringToFileMap-0", "gs://path/to/stringTofile1", DefaultPathBuilder.get("path/to/stringTofile1"), workingDisk))
        jesInputs should contain(JesFileInput(
          "stringToFileMap-1", "gs://path/to/stringTofile2", DefaultPathBuilder.get("path/to/stringTofile2"), workingDisk))
        jesInputs should contain(JesFileInput(
          "fileToStringMap-0", "gs://path/to/fileToString1", DefaultPathBuilder.get("path/to/fileToString1"), workingDisk))
        jesInputs should contain(JesFileInput(
          "fileToStringMap-1", "gs://path/to/fileToString2", DefaultPathBuilder.get("path/to/fileToString2"), workingDisk))
        jesInputs should contain(JesFileInput(
          "fileToFileMap-0", "gs://path/to/fileToFile1Key", DefaultPathBuilder.get("path/to/fileToFile1Key"), workingDisk))
        jesInputs should contain(JesFileInput(
          "fileToFileMap-1", "gs://path/to/fileToFile1Value", DefaultPathBuilder.get("path/to/fileToFile1Value"), workingDisk))
        jesInputs should contain(JesFileInput(
          "fileToFileMap-2", "gs://path/to/fileToFile2Key", DefaultPathBuilder.get("path/to/fileToFile2Key"), workingDisk))
        jesInputs should contain(JesFileInput(
          "fileToFileMap-3", "gs://path/to/fileToFile2Value", DefaultPathBuilder.get("path/to/fileToFile2Value"), workingDisk))

      case Left(badness) => fail(badness.toList.mkString(", "))
    }
  }

  def makeJesActorRef(sampleWdl: SampleWdl, callName: LocallyQualifiedName, inputs: Map[FullyQualifiedName, WomValue],
                      functions: JesExpressionFunctions = TestableJesExpressionFunctions):
  TestActorRef[TestableJesJobExecutionActor] = {
    val womWorkflow = WdlNamespaceWithWorkflow.load(sampleWdl.asWorkflowSources(DockerAndDiskRuntime).workflowSource,
      Seq.empty[ImportResolver]).get.workflow.toWomWorkflowDefinition.getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    dockerAndDiskWdlNamespace.toWomExecutable(Option(inputs.toJson.compactPrint)) match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap({case (port, v) => v.select[WomValue] map { port -> _ }})
        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty
        )

        val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.find(_.localName == callName).get
        val key = BackendJobDescriptorKey(call, None, 1)
        val runtimeAttributes = makeRuntimeAttributes(call)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, fqnWdlMapToDeclarationMap(inputs), NoDocker, Map.empty)

        val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration, functions))
        TestActorRef[TestableJesJobExecutionActor](props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")
      case Left(badness) => fail(badness.toList.mkString(", "))
    }
  }

  it should "generate correct JesOutputs" taggedAs PostWomTest ignore {
    val inputs = Map(
      "in" -> WomSingleFile("gs://blah/b/c.txt")
    )
    val jesBackend = makeJesActorRef(SampleWdl.FilePassingWorkflow, "a", inputs).underlyingActor
    val jobDescriptor = jesBackend.jobDescriptor
    val workflowId = jesBackend.workflowId
    val jesInputs = jesBackend.generateJesInputs(jobDescriptor)
    jesInputs should have size 1
    jesInputs should contain(JesFileInput("in-0", "gs://blah/b/c.txt", DefaultPathBuilder.get("blah/b/c.txt"), workingDisk))
    val jesOutputs = jesBackend.generateJesOutputs(jobDescriptor)
    jesOutputs should have size 1
    jesOutputs should contain(JesFileOutput("out",
      s"gs://my-cromwell-workflows-bucket/file_passing/$workflowId/call-a/out", DefaultPathBuilder.get("out"), workingDisk))
  }

  it should "generate correct JesInputs when a command line contains a write_lines call in it" taggedAs PostWomTest ignore {
    val inputs = Map(
      "strs" -> WomArray(WomArrayType(WomStringType), Seq("A", "B", "C").map(WomString))
    )

    class TestJesExpressionFunctions extends JesExpressionFunctions(TestableStandardExpressionFunctionsParams) {
      override def writeFile(path: String, content: String): Future[WomSingleFile] = {
        Future.fromTry(Success(WomSingleFile(s"gs://some/path/file.txt")))
      }
    }

    val functions = new TestJesExpressionFunctions
    val jesBackend = makeJesActorRef(SampleWdl.ArrayIO, "serialize", inputs, functions).underlyingActor
    val jobDescriptor = jesBackend.jobDescriptor
    val jesInputs = jesBackend.generateJesInputs(jobDescriptor)
    jesInputs should have size 1
    jesInputs should contain(JesFileInput(
      "c6fd5c91-0", "gs://some/path/file.txt", DefaultPathBuilder.get("some/path/file.txt"), workingDisk))
    val jesOutputs = jesBackend.generateJesOutputs(jobDescriptor)
    jesOutputs should have size 0
  }

  it should "generate correct JesFileInputs from a WdlArray" taggedAs PostWomTest ignore {
    val inputs: Map[String, WomValue] = Map(
      "fileArray" ->
        WomArray(WomArrayType(WomSingleFileType), Seq(WomSingleFile("gs://path/to/file1"), WomSingleFile("gs://path/to/file2")))
    )

    val womWorkflow = dockerAndDiskWdlNamespace.workflow.toWomWorkflowDefinition.getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    dockerAndDiskWdlNamespace.toWomExecutable(Option(inputs.toJson.compactPrint)) match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap({case (port, v) => v.select[WomValue] map { port -> _ }})
        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty
        )

        val job: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
        val runtimeAttributes = makeRuntimeAttributes(job)
        val key = BackendJobDescriptorKey(job, None, 1)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, fqnWdlMapToDeclarationMap(inputs), NoDocker, Map.empty)

        val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
        val testActorRef = TestActorRef[TestableJesJobExecutionActor](
          props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

        val jesInputs = testActorRef.underlyingActor.generateJesInputs(jobDescriptor)
        jesInputs should have size 2
        jesInputs should contain(JesFileInput("fileArray-0", "gs://path/to/file1", DefaultPathBuilder.get("path/to/file1"), workingDisk))
        jesInputs should contain(JesFileInput("fileArray-1", "gs://path/to/file2", DefaultPathBuilder.get("path/to/file2"), workingDisk))
      case Left(badness) => fail(badness.toList.mkString(", "))
    }
  }

  it should "generate correct JesFileInputs from a WdlFile" taggedAs PostWomTest ignore {
    val inputs: Map[String, WomValue] = Map(
      "file1" -> WomSingleFile("gs://path/to/file1"),
      "file2" -> WomSingleFile("gs://path/to/file2")
    )

    val womWorkflow = dockerAndDiskWdlNamespace.workflow.toWomWorkflowDefinition.getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    dockerAndDiskWdlNamespace.toWomExecutable(Option(inputs.toJson.compactPrint)) match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap({case (port, v) => v.select[WomValue] map { port -> _ }})
        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty
        )

        val job: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
        val runtimeAttributes = makeRuntimeAttributes(job)
        val key = BackendJobDescriptorKey(job, None, 1)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, fqnWdlMapToDeclarationMap(inputs), NoDocker, Map.empty)

        val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
        val testActorRef = TestActorRef[TestableJesJobExecutionActor](
          props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

        val jesInputs = testActorRef.underlyingActor.generateJesInputs(jobDescriptor)
        jesInputs should have size 2
        jesInputs should contain(JesFileInput("file1-0", "gs://path/to/file1", DefaultPathBuilder.get("path/to/file1"), workingDisk))
        jesInputs should contain(JesFileInput("file2-0", "gs://path/to/file2", DefaultPathBuilder.get("path/to/file2"), workingDisk))

      case Left(badness) => fail(badness.toList.mkString(", "))
    }
  }

  it should "convert local Paths back to corresponding GCS paths in JesOutputs" in {
    val jesOutputs = Set(
      JesFileOutput("/cromwell_root/path/to/file1", "gs://path/to/file1",
        DefaultPathBuilder.get("/cromwell_root/path/to/file1"), workingDisk),
      JesFileOutput("/cromwell_root/path/to/file2", "gs://path/to/file2",
        DefaultPathBuilder.get("/cromwell_root/path/to/file2"), workingDisk),
      JesFileOutput("/cromwell_root/path/to/file3", "gs://path/to/file3",
        DefaultPathBuilder.get("/cromwell_root/path/to/file3"), workingDisk),
      JesFileOutput("/cromwell_root/path/to/file4", "gs://path/to/file4",
        DefaultPathBuilder.get("/cromwell_root/path/to/file4"), workingDisk),
      JesFileOutput("/cromwell_root/path/to/file5", "gs://path/to/file5",
        DefaultPathBuilder.get("/cromwell_root/path/to/file5"), workingDisk)
    )
    val outputValues = Seq(
      WomSingleFile("/cromwell_root/path/to/file1"),
      WomArray(WomArrayType(WomSingleFileType), Seq(
        WomSingleFile("/cromwell_root/path/to/file2"), WomSingleFile("/cromwell_root/path/to/file3"))),
      WomMap(WomMapType(WomSingleFileType, WomSingleFileType), Map(
        WomSingleFile("/cromwell_root/path/to/file4") -> WomSingleFile("/cromwell_root/path/to/file5")
      ))
    )

    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      WdlNamespaceWithWorkflow.load(SampleWdl.EmptyString.asWorkflowSources(DockerAndDiskRuntime).workflowSource,
        Seq.empty[ImportResolver]).get.workflow.toWomWorkflowDefinition.getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      NoOptions,
      Labels.empty
    )

    val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
    val key = BackendJobDescriptorKey(call, None, 1)
    val runtimeAttributes = makeRuntimeAttributes(call)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, Map.empty)

    val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
    val testActorRef = TestActorRef[TestableJesJobExecutionActor](
      props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

    def wdlValueToGcsPath(jesOutputs: Set[JesFileOutput])(womValue: WomValue): WomValue = {
      WomFileMapper.mapWomFiles(testActorRef.underlyingActor.womFileToGcsPath(jesOutputs))(womValue).get
    }

    val result = outputValues map wdlValueToGcsPath(jesOutputs)
    result should have size 3
    result should contain(WomSingleFile("gs://path/to/file1"))
    result should contain(WomArray(WomArrayType(WomSingleFileType),
      Seq(WomSingleFile("gs://path/to/file2"), WomSingleFile("gs://path/to/file3")))
    )
    result should contain(WomMap(WomMapType(WomSingleFileType, WomSingleFileType),
      Map(WomSingleFile("gs://path/to/file4") -> WomSingleFile("gs://path/to/file5")))
    )
  }

  it should "create a JesFileInput for the monitoring script, when specified" in {
    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      WdlNamespaceWithWorkflow.load(SampleWdl.EmptyString.asWorkflowSources(DockerAndDiskRuntime).workflowSource,
        Seq.empty[ImportResolver]).get.workflow.toWomWorkflowDefinition.getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      WorkflowOptions.fromJsonString("""{"monitoring_script": "gs://path/to/script"}""").get,
      Labels.empty
    )

    val job: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
    val runtimeAttributes = makeRuntimeAttributes(job)
    val key = BackendJobDescriptorKey(job, None, 1)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, Map.empty)

    val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
    val testActorRef = TestActorRef[TestableJesJobExecutionActor](
      props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

    testActorRef.underlyingActor.monitoringScript shouldBe
      Some(JesFileInput("monitoring-in", "gs://path/to/script", DefaultPathBuilder.get("monitoring.sh"), workingDisk))
  }

  it should "not create a JesFileInput for the monitoring script, when not specified" in {
    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      WdlNamespaceWithWorkflow.load(SampleWdl.EmptyString.asWorkflowSources(DockerAndDiskRuntime).workflowSource,
        Seq.empty[ImportResolver]).get.workflow.toWomWorkflowDefinition.getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      NoOptions,
      Labels.empty
    )

    val job: CommandCallNode = workflowDescriptor.callable.graph.nodes.collectFirst({case t: CommandCallNode => t}).get
    val key = BackendJobDescriptorKey(job, None, 1)
    val runtimeAttributes = makeRuntimeAttributes(job)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, Map.empty)

    val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
    val testActorRef = TestActorRef[TestableJesJobExecutionActor](
      props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

    testActorRef.underlyingActor.monitoringScript shouldBe None
  }

  it should "return JES log paths for non-scattered call" in {
    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId(UUID.fromString("e6236763-c518-41d0-9688-432549a8bf7c")),
      WdlNamespaceWithWorkflow.load(
        SampleWdl.HelloWorld.asWorkflowSources(""" runtime {docker: "ubuntu:latest"} """).workflowSource,
        Seq.empty[ImportResolver]).get.workflow.toWomWorkflowDefinition.getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      WorkflowOptions.fromJsonString(""" {"jes_gcs_root": "gs://path/to/gcs_root"} """).get,
      Labels.empty
    )

    val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.find(_.localName == "hello").get
    val key = BackendJobDescriptorKey(call, None, 1)
    val runtimeAttributes = makeRuntimeAttributes(call)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, Map.empty)

    val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
    val testActorRef = TestActorRef[TestableJesJobExecutionActor](
      props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

    val jesBackend = testActorRef.underlyingActor

    jesBackend.jesCallPaths.stdout should be(a[GcsPath])
    jesBackend.jesCallPaths.stdout.pathAsString shouldBe
      "gs://path/to/gcs_root/wf_hello/e6236763-c518-41d0-9688-432549a8bf7c/call-hello/hello-stdout.log"
    jesBackend.jesCallPaths.stderr should be(a[GcsPath])
    jesBackend.jesCallPaths.stderr.pathAsString shouldBe
      "gs://path/to/gcs_root/wf_hello/e6236763-c518-41d0-9688-432549a8bf7c/call-hello/hello-stderr.log"
    jesBackend.jesCallPaths.jesLogPath should be(a[GcsPath])
    jesBackend.jesCallPaths.jesLogPath.pathAsString shouldBe
      "gs://path/to/gcs_root/wf_hello/e6236763-c518-41d0-9688-432549a8bf7c/call-hello/hello.log"
  }

  it should "return JES log paths for scattered call" taggedAs PostWomTest ignore {
    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId(UUID.fromString("e6236763-c518-41d0-9688-432549a8bf7d")),
      WdlNamespaceWithWorkflow.load(
        new SampleWdl.ScatterWdl().asWorkflowSources(""" runtime {docker: "ubuntu:latest"} """).workflowSource,
        Seq.empty[ImportResolver]).get.workflow.toWomWorkflowDefinition.getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      WorkflowOptions.fromJsonString(""" {"jes_gcs_root": "gs://path/to/gcs_root"} """).get,
      Labels.empty
    )

    val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.find(_.localName == "B").get
    val key = BackendJobDescriptorKey(call, Option(2), 1)
    val runtimeAttributes = makeRuntimeAttributes(call)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, Map.empty)

    val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
    val testActorRef = TestActorRef[TestableJesJobExecutionActor](
      props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

    val jesBackend = testActorRef.underlyingActor

    jesBackend.jesCallPaths.stdout should be(a[GcsPath])
    jesBackend.jesCallPaths.stdout.pathAsString shouldBe
      "gs://path/to/gcs_root/w/e6236763-c518-41d0-9688-432549a8bf7d/call-B/shard-2/B-2-stdout.log"
    jesBackend.jesCallPaths.stderr should be(a[GcsPath])
    jesBackend.jesCallPaths.stderr.pathAsString shouldBe
      "gs://path/to/gcs_root/w/e6236763-c518-41d0-9688-432549a8bf7d/call-B/shard-2/B-2-stderr.log"
    jesBackend.jesCallPaths.jesLogPath should be(a[GcsPath])
    jesBackend.jesCallPaths.jesLogPath.pathAsString shouldBe
      "gs://path/to/gcs_root/w/e6236763-c518-41d0-9688-432549a8bf7d/call-B/shard-2/B-2.log"
  }

  it should "return preemptible = true only in the correct cases" in {
    def attempt(max: Int, attempt: Int): JesAsyncBackendJobExecutionActor = {
      buildPreemptibleTestActorRef(attempt, max).underlyingActor
    }
    def attempt1(max: Int) = attempt(max, 1)
    def attempt2(max: Int) = attempt(max, 2)

    val descriptorWithMax0AndKey1 = attempt1(max = 0)
    descriptorWithMax0AndKey1.preemptible shouldBe false

    val descriptorWithMax1AndKey1 = attempt1(max = 1)
    descriptorWithMax1AndKey1.preemptible shouldBe true

    val descriptorWithMax2AndKey1 = attempt1(max = 2)
    descriptorWithMax2AndKey1.preemptible shouldBe true

    val descriptorWithMax1AndKey2 = attempt2(max = 1)
    descriptorWithMax1AndKey2.preemptible shouldBe false

    val descriptorWithMax2AndKey2 = attempt2(max = 2)
    descriptorWithMax2AndKey2.preemptible shouldBe true
  }

  private def makeRuntimeAttributes(job: CommandCallNode) = {
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(job.callable.runtimeAttributes, null, Map.empty)
    RuntimeAttributeDefinition.addDefaultsToAttributes(
      runtimeAttributesBuilder.definitions.toSet, NoOptions)(evaluatedAttributes.getOrElse(fail("Failed to evaluate runtime attributes")))
  }
}
