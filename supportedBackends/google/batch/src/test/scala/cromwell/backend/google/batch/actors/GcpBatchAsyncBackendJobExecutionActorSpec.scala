package cromwell.backend.google.batch
package actors

import _root_.wdl.draft2.model._
import akka.actor.{ActorRef, Props}
import akka.testkit.{ImplicitSender, TestActorRef, TestDuration, TestProbe}
import cats.data.NonEmptyList
import cloud.nio.impl.drs.DrsCloudNioFileProvider.DrsReadInterpreter
import cloud.nio.impl.drs.{DrsCloudNioFileSystemProvider, GoogleOauthDrsCredentials}
import com.google.cloud.NoCredentials
import com.google.cloud.batch.v1.{CreateJobRequest, DeleteJobRequest, GetJobRequest, JobName}
import com.typesafe.config.{Config, ConfigFactory}
import common.collections.EnhancedCollections._
import common.mock.MockSugar
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend._
import cromwell.backend.async.AsyncBackendJobExecutionActor.{Execute, ExecutionMode}
import cromwell.backend.async.{ExecutionHandle, FailedNonRetryableExecutionHandle, FailedRetryableExecutionHandle}
import cromwell.backend.google.batch.actors.GcpBatchAsyncBackendJobExecutionActor.GcpBatchPendingExecutionHandle
import cromwell.backend.google.batch.api.GcpBatchRequestFactory
import cromwell.backend.google.batch.io.{DiskType, GcpBatchWorkingDisk}
import cromwell.backend.google.batch.models._
import cromwell.backend.google.batch.runnable.RunnableUtils.MountPoint
import cromwell.backend.google.batch.util.BatchExpressionFunctions
import cromwell.backend.io.JobPathsSpecHelper._
import cromwell.backend.standard.{
  DefaultStandardAsyncExecutionActorParams,
  StandardAsyncExecutionActorParams,
  StandardAsyncJob,
  StandardExpressionFunctionsParams
}
import cromwell.core._
import cromwell.core.callcaching.NoDocker
import cromwell.core.labels.Labels
import cromwell.core.logging.JobLogger
import cromwell.core.path.{DefaultPathBuilder, PathBuilder}
import cromwell.filesystems.drs.DrsPathBuilder
import cromwell.filesystems.gcs.{GcsPath, GcsPathBuilder, MockGcsPathBuilder}
import cromwell.services.instrumentation.InstrumentationService.InstrumentationServiceMessage
import cromwell.services.instrumentation.{CromwellBucket, CromwellIncrement}
import cromwell.services.keyvalue.InMemoryKvServiceActor
import cromwell.services.keyvalue.KeyValueServiceActor.{KvJobKey, KvPair, ScopedKey}
import cromwell.services.metadata.CallMetadataKeys
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metrics.bard.BardEventing.BardEventRequest
import cromwell.services.metrics.bard.model.TaskSummaryEvent
import cromwell.util.JsonFormatting.WomValueJsonFormatter._
import cromwell.util.SampleWdl
import org.mockito.Mockito._
import org.scalatest._
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.slf4j.Logger
import spray.json._
import wdl.transforms.draft2.wdlom2wom.WdlDraft2WomExecutableMakers._
import wdl.transforms.draft2.wdlom2wom._
import wom.WomFileMapper
import wom.expression.NoIoFunctionSet
import wom.graph.CommandCallNode
import wom.transforms.WomExecutableMaker.ops._
import wom.transforms.WomWorkflowDefinitionMaker.ops._
import wom.types.{WomArrayType, WomMapType, WomSingleFileType, WomStringType}
import wom.values._

import java.nio.file.Paths
import java.time.OffsetDateTime
import java.time.temporal.ChronoUnit
import java.util.UUID
import scala.concurrent.duration._
import scala.concurrent.{Await, ExecutionContext, Future, Promise}
import scala.language.postfixOps
import scala.util.Success

class GcpBatchAsyncBackendJobExecutionActorSpec
    extends TestKitSuite
    with AnyFlatSpecLike
    with Matchers
    with ImplicitSender
    with BackendSpec
    with BeforeAndAfter
    with MockSugar
    with DefaultJsonProtocol {

  val mockPathBuilder: GcsPathBuilder = MockGcsPathBuilder.instance
  import MockGcsPathBuilder._
  var kvService: ActorRef = system.actorOf(Props(new InMemoryKvServiceActor), "kvService")

  private def gcsPath(str: String) = mockPathBuilder.build(str).getOrElse(fail(s"Invalid gcs path: $str"))

  import cromwell.backend.google.batch.models.GcpBatchTestConfig._

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

  private val NoOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue]))

  private lazy val TestableCallContext =
    CallContext(mockPathBuilder.build("gs://root").get, DummyStandardPaths, isDocker = false)

  private lazy val TestableStandardExpressionFunctionsParams: StandardExpressionFunctionsParams =
    new StandardExpressionFunctionsParams {
      override lazy val pathBuilders: List[PathBuilder] = List(mockPathBuilder)
      override lazy val callContext: CallContext = TestableCallContext
      override val ioActorProxy: ActorRef = simpleIoActor
      override val executionContext: ExecutionContext = system.dispatcher
    }

  lazy val TestableGcpBatchExpressionFunctions: BatchExpressionFunctions =
    new BatchExpressionFunctions(TestableStandardExpressionFunctionsParams)

  private def buildInitializationData(jobDescriptor: BackendJobDescriptor, configuration: GcpBatchConfiguration) = {
    val pathBuilders =
      Await.result(configuration.configurationDescriptor.pathBuilders(WorkflowOptions.empty), 5.seconds)
    val workflowPaths = GcpBatchWorkflowPaths(
      jobDescriptor.workflowDescriptor,
      NoCredentials.getInstance(),
      NoCredentials.getInstance(),
      configuration,
      pathBuilders,
      GcpBatchInitializationActor.defaultStandardStreamNameToFileNameMetadataMapper
    )
    val runtimeAttributesBuilder = GcpBatchRuntimeAttributes.runtimeAttributesBuilder(configuration)

    val requestFactory: GcpBatchRequestFactory = new GcpBatchRequestFactory {
      override def submitRequest(data: GcpBatchRequest, jobLogger: JobLogger): CreateJobRequest = null

      override def queryRequest(jobName: JobName): GetJobRequest = null

      override def abortRequest(jobName: JobName): DeleteJobRequest = null
    }
    GcpBackendInitializationData(workflowPaths,
                                 runtimeAttributesBuilder,
                                 configuration,
                                 null,
                                 None,
                                 None,
                                 None,
                                 requestFactory
    )
  }

  class TestableGcpBatchJobExecutionActor(params: StandardAsyncExecutionActorParams,
                                          functions: BatchExpressionFunctions
  ) extends GcpBatchAsyncBackendJobExecutionActor(params) {

    def this(jobDescriptor: BackendJobDescriptor,
             promise: Promise[BackendJobExecutionResponse],
             batchConfiguration: GcpBatchConfiguration,
             functions: BatchExpressionFunctions = TestableGcpBatchExpressionFunctions,
             batchSingletonActor: ActorRef = emptyActor,
             ioActor: ActorRef = mockIoActor,
             serviceRegistryActor: ActorRef = kvService
    ) =
      this(
        DefaultStandardAsyncExecutionActorParams(
          jobIdKey = GcpBatchAsyncBackendJobExecutionActor.GcpBatchOperationIdKey,
          serviceRegistryActor = serviceRegistryActor,
          ioActor = ioActor,
          jobDescriptor = jobDescriptor,
          configurationDescriptor = batchConfiguration.configurationDescriptor,
          backendInitializationDataOption = Option(buildInitializationData(jobDescriptor, batchConfiguration)),
          backendSingletonActorOption = Option(batchSingletonActor),
          completionPromise = promise,
          minimumRuntimeSettings = MinimumRuntimeSettings(),
          groupMetricsActor = emptyActor
        ),
        functions
      )

    override lazy val jobLogger: JobLogger = new JobLogger(
      loggerName = "TestLogger",
      workflowIdForLogging = workflowId.toPossiblyNotRoot,
      rootWorkflowIdForLogging = workflowId.toRoot,
      jobTag = jobTag,
      akkaLogger = Option(log)
    ) {
      override def tag: String = s"$name [UUID(${workflowId.shortString})$jobTag]"
      override val slf4jLoggers: Set[Logger] = Set.empty
    }

    override lazy val backendEngineFunctions: BatchExpressionFunctions = functions

    override val pollingResultMonitorActor: Option[ActorRef] = Some(
      context.actorOf(
        BatchPollResultMonitorActor.props(serviceRegistryActor,
                                          workflowDescriptor,
                                          jobDescriptor,
                                          validatedRuntimeAttributes,
                                          platform,
                                          jobLogger
        )
      )
    )
  }

  private val runtimeAttributesBuilder = GcpBatchRuntimeAttributes.runtimeAttributesBuilder(gcpBatchConfiguration)
  private val workingDisk = GcpBatchWorkingDisk(DiskType.SSD, 200)

  val DockerAndDiskRuntime: String =
    """
      |runtime {
      |  docker: "ubuntu:latest"
      |  disks: "local-disk 200 SSD"
      |}
    """.stripMargin

  private def buildPreemptibleJobDescriptor(preemptible: Int,
                                            previousPreemptions: Int,
                                            previousUnexpectedRetries: Int,
                                            previousTransientRetries: Int,
                                            failedRetriesCountOpt: Option[Int] = None
  ): BackendJobDescriptor = {
    val attempt = previousPreemptions + previousUnexpectedRetries + 1
    val wdlNamespace = WdlNamespaceWithWorkflow
      .load(YoSup.replace("[PREEMPTIBLE]", s"preemptible: $preemptible"), Seq.empty[Draft2ImportResolver])
      .get
    val womDefinition = wdlNamespace.workflow
      .toWomWorkflowDefinition(isASubworkflow = false)
      .getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))

    wdlNamespace.toWomExecutable(Option(Inputs.toJson.compactPrint), NoIoFunctionSet, strictValidation = true) match {
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
          Labels.empty,
          HogGroup("foo"),
          List.empty,
          None
        )

        val job = workflowDescriptor.callable.taskCallNodes.head
        val key = BackendJobDescriptorKey(job, None, attempt)
        val runtimeAttributes = makeRuntimeAttributes(job)
        val prefetchedKvEntries = Map(
          GcpBatchBackendLifecycleActorFactory.preemptionCountKey -> KvPair(
            ScopedKey(workflowDescriptor.id, KvJobKey(key), GcpBatchBackendLifecycleActorFactory.preemptionCountKey),
            previousPreemptions.toString
          ),
          GcpBatchBackendLifecycleActorFactory.unexpectedRetryCountKey -> KvPair(
            ScopedKey(workflowDescriptor.id,
                      KvJobKey(key),
                      GcpBatchBackendLifecycleActorFactory.unexpectedRetryCountKey
            ),
            previousUnexpectedRetries.toString
          ),
          GcpBatchBackendLifecycleActorFactory.transientRetryCountKey -> KvPair(
            ScopedKey(workflowDescriptor.id,
                      KvJobKey(key),
                      GcpBatchBackendLifecycleActorFactory.transientRetryCountKey
            ),
            previousTransientRetries.toString
          )
        )
        val prefetchedKvEntriesUpd = if (failedRetriesCountOpt.isEmpty) {
          prefetchedKvEntries
        } else {
          prefetchedKvEntries + (BackendLifecycleActorFactory.FailedRetryCountKey -> KvPair(
            ScopedKey(workflowDescriptor.id, KvJobKey(key), BackendLifecycleActorFactory.FailedRetryCountKey),
            failedRetriesCountOpt.get.toString
          ))
        }
        BackendJobDescriptor(workflowDescriptor,
                             key,
                             runtimeAttributes,
                             fqnWdlMapToDeclarationMap(Inputs),
                             NoDocker,
                             None,
                             prefetchedKvEntriesUpd
        )
      case Left(badtimes) => fail(badtimes.toList.mkString(", "))
    }
  }

  private def executionActor(jobDescriptor: BackendJobDescriptor,
                             promise: Promise[BackendJobExecutionResponse],
                             batchSingletonActor: ActorRef,
                             shouldBePreemptible: Boolean,
                             serviceRegistryActor: ActorRef,
                             referenceInputFilesOpt: Option[Set[GcpBatchInput]]
  ): ActorRef = {

    val job = generateStandardAsyncJob
    val run = Run(job)
    val handle = new GcpBatchPendingExecutionHandle(jobDescriptor, run.job, Option(run), None)

    class ExecuteOrRecoverActor
        extends TestableGcpBatchJobExecutionActor(jobDescriptor,
                                                  promise,
                                                  gcpBatchConfiguration,
                                                  batchSingletonActor = batchSingletonActor,
                                                  serviceRegistryActor = serviceRegistryActor
        ) {
      override def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
        sendIncrementMetricsForReferenceFiles(referenceInputFilesOpt)

        if (preemptible == shouldBePreemptible) Future.successful(handle)
        else Future.failed(new Exception(s"Test expected preemptible to be $shouldBePreemptible but got $preemptible"))
      }
    }

    system.actorOf(Props(new ExecuteOrRecoverActor), "ExecuteOrRecoverActor-" + UUID.randomUUID)
  }

  def buildPreemptibleTestActorRef(attempt: Int,
                                   preemptible: Int,
                                   previousTransientRetriesCount: Int = 0,
                                   failedRetriesCountOpt: Option[Int] = None
  ): TestActorRef[TestableGcpBatchJobExecutionActor] = {
    // For this test we say that all previous attempts were preempted:
    val jobDescriptor = buildPreemptibleJobDescriptor(preemptible,
                                                      attempt - 1,
                                                      previousUnexpectedRetries = 0,
                                                      previousTransientRetries = previousTransientRetriesCount,
                                                      failedRetriesCountOpt = failedRetriesCountOpt
    )
    val props = Props(
      new TestableGcpBatchJobExecutionActor(jobDescriptor,
                                            Promise(),
                                            gcpBatchConfiguration,
                                            TestableGcpBatchExpressionFunctions,
                                            emptyActor,
                                            failIoActor
      )
    )
    TestActorRef(props, s"TestableGcpBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")
  }

  behavior of "GcpBatchAsyncBackendJobExecutionActor"

  private val timeout = 25 seconds

  it should "group files by bucket" in {

    def makeInput(bucket: String, name: String): GcpBatchFileInput = {
      val mockCloudPath = mock[cromwell.core.path.Path]
      when(mockCloudPath.pathAsString) thenReturn s"gs://$bucket/$name"

      GcpBatchFileInput(
        name = name,
        cloudPath = mockCloudPath,
        relativeHostPath = DefaultPathBuilder.build(Paths.get(s"$bucket/$name")),
        mount = null
      )
    }

    val inputs: List[GcpBatchFileInput] = List(
      ("foo", "file1"),
      ("foo", "file2"),
      ("bar", "file1"),
      ("bar", "file2"),
      ("bar", "file3"),
      ("baz", "file1")
    ) map (makeInput _).tupled.apply

    val expected =
      Map("foo" -> (NonEmptyList.of(0, 1) map inputs.apply)) ++
        Map("bar" -> (NonEmptyList.of(2, 3, 4) map inputs.apply)) ++
        Map("baz" -> NonEmptyList.of(inputs(5)))

    GcpBatchAsyncBackendJobExecutionActor.groupParametersByGcsBucket(inputs) shouldEqual expected
  }

  it should "generate a CSV manifest for DRS inputs, ignoring non-DRS inputs" in {
    def makeDrsPathBuilder: DrsPathBuilder = {
      val drsResolverConfig: Config = ConfigFactory.parseString(
        """resolver {
          |   url = "http://drshub-url"
          |}
          |""".stripMargin
      )

      val fakeCredentials = NoCredentials.getInstance

      val drsReadInterpreter: DrsReadInterpreter = (_, _) =>
        throw new UnsupportedOperationException(
          "GcpBatchAsyncBackendJobExecutionActorSpec doesn't need to use drs read interpreter."
        )

      DrsPathBuilder(
        new DrsCloudNioFileSystemProvider(drsResolverConfig,
                                          GoogleOauthDrsCredentials(fakeCredentials, 1.minutes),
                                          drsReadInterpreter
        ),
        None
      )
    }

    val mount = GcpBatchWorkingDisk(DiskType.LOCAL, 1)

    def makeDrsInput(name: String, drsUri: String, containerPath: String): GcpBatchFileInput = {
      val drsPath = makeDrsPathBuilder.build(drsUri).get
      val containerRelativePath = DefaultPathBuilder.get(containerPath)
      GcpBatchFileInput(name, drsPath, containerRelativePath, mount)
    }

    val nonDrsInput: GcpBatchFileInput = GcpBatchFileInput("nnn",
                                                           DefaultPathBuilder.get("/local/nnn.bai"),
                                                           DefaultPathBuilder.get("/path/to/nnn.bai"),
                                                           mount
    )

    val inputs = List(
      makeDrsInput("aaa", "drs://drs.example.org/aaa", "path/to/aaa.bai"),
      nonDrsInput,
      makeDrsInput("bbb", "drs://drs.example.org/bbb", "path/to/bbb.bai")
    )

    GcpBatchAsyncBackendJobExecutionActor.generateDrsLocalizerManifest(inputs) shouldEqual
      s"drs://drs.example.org/aaa,$MountPoint/path/to/aaa.bai\r\ndrs://drs.example.org/bbb,$MountPoint/path/to/bbb.bai\r\n"
  }

  it should "send proper value for \"number of reference files used gauge\" metric, or don't send anything if reference disks feature is disabled" in {

    val expectedInput1 = GcpBatchFileInput(name = "testfile1",
                                           relativeHostPath =
                                             DefaultPathBuilder.build(Paths.get(s"test/reference/path/file1")),
                                           mount = null,
                                           cloudPath = null
    )
    val expectedInput2 = GcpBatchFileInput(name = "testfile2",
                                           relativeHostPath =
                                             DefaultPathBuilder.build(Paths.get(s"test/reference/path/file2")),
                                           mount = null,
                                           cloudPath = null
    )
    val expectedReferenceInputFiles = Set[GcpBatchInput](expectedInput1, expectedInput2)

    val expectedMsg1 = InstrumentationServiceMessage(
      CromwellIncrement(
        CromwellBucket(List.empty, NonEmptyList.of("referencefiles", expectedInput1.relativeHostPath.pathAsString))
      )
    )
    val expectedMsg2 = InstrumentationServiceMessage(
      CromwellIncrement(
        CromwellBucket(List.empty, NonEmptyList.of("referencefiles", expectedInput2.relativeHostPath.pathAsString))
      )
    )

    val jobDescriptor = buildPreemptibleJobDescriptor(0, 0, 0, 0)
    val serviceRegistryProbe = TestProbe()

    val backend1 = executionActor(
      jobDescriptor,
      Promise[BackendJobExecutionResponse](),
      TestProbe().ref,
      shouldBePreemptible = false,
      serviceRegistryActor = serviceRegistryProbe.ref,
      referenceInputFilesOpt = Option(expectedReferenceInputFiles)
    )
    backend1 ! Execute
    serviceRegistryProbe.expectMsgAllOf(expectedMsg1, expectedMsg2)

    val backend2 = executionActor(
      jobDescriptor,
      Promise[BackendJobExecutionResponse](),
      TestProbe().ref,
      shouldBePreemptible = false,
      serviceRegistryActor = serviceRegistryProbe.ref,
      referenceInputFilesOpt = None
    )
    backend2 ! Execute
    serviceRegistryProbe.expectNoMessage(timeout)
  }

  it should "not restart 2 of 1 unexpected shutdowns without another preemptible VM" in {

    val actorRef = buildPreemptibleTestActorRef(2, 1, 0)
    val batchBackend = actorRef.underlyingActor
    val runId = generateStandardAsyncJob
    val handle = new GcpBatchPendingExecutionHandle(null, runId, None, None)

    val failedStatus = RunStatus.Failed(
      GcpBatchExitCode.Success,
      Seq.empty
    )
    val executionResult = batchBackend.handleExecutionResult(failedStatus, handle)
    val result = Await.result(executionResult, timeout)
    result.isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    val failedHandle = result.asInstanceOf[FailedNonRetryableExecutionHandle]
    failedHandle.returnCode shouldBe None
  }

  it should "retry transient failures when appropriate" in {
    val actorRef = buildPreemptibleTestActorRef(attempt = 2, preemptible = 0, previousTransientRetriesCount = 1)
    val batchBackend = actorRef.underlyingActor
    val runId = generateStandardAsyncJob
    val handle = new GcpBatchPendingExecutionHandle(null, runId, None, None)

    val failedStatus = RunStatus.Failed(
      GcpBatchExitCode.VMRecreatedDuringExecution,
      Seq.empty
    )
    val executionResult = batchBackend.handleExecutionResult(failedStatus, handle)
    val result = Await.result(executionResult, timeout)
    result.isInstanceOf[FailedRetryableExecutionHandle] shouldBe true
    val failedHandle = result.asInstanceOf[FailedRetryableExecutionHandle]
    failedHandle.returnCode shouldBe None
    failedHandle.kvPairsToSave.map { pairs =>
      pairs.exists(p => p.key.key == GcpBatchBackendLifecycleActorFactory.preemptionCountKey && p.value == "0")
      pairs.exists(p => p.key.key == GcpBatchBackendLifecycleActorFactory.unexpectedRetryCountKey && p.value == "0")
      pairs.exists(p => p.key.key == GcpBatchBackendLifecycleActorFactory.transientRetryCountKey && p.value == "2")
    }
  }

  it should "not retry transient failures after 10 attempts" in {
    val actorRef = buildPreemptibleTestActorRef(attempt = 11, preemptible = 0, previousTransientRetriesCount = 10)
    val batchBackend = actorRef.underlyingActor
    val runId = generateStandardAsyncJob
    val handle = new GcpBatchPendingExecutionHandle(null, runId, None, None)

    val failedStatus = RunStatus.Failed(
      GcpBatchExitCode.VMRecreatedDuringExecution,
      Seq.empty
    )
    val executionResult = batchBackend.handleExecutionResult(failedStatus, handle)
    val result = Await.result(executionResult, timeout)
    result.isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    val failedHandle = result.asInstanceOf[FailedNonRetryableExecutionHandle]
    failedHandle.returnCode shouldBe None
    failedHandle.kvPairsToSave.map { pairs =>
      pairs.exists(p => p.key.key == GcpBatchBackendLifecycleActorFactory.preemptionCountKey && p.value == "0")
      pairs.exists(p => p.key.key == GcpBatchBackendLifecycleActorFactory.unexpectedRetryCountKey && p.value == "1")
      pairs.exists(p => p.key.key == GcpBatchBackendLifecycleActorFactory.transientRetryCountKey && p.value == "10")
    }
  }

  it should "choose transient retry over preemptible retry when task has not started" in {
    val actorRef = buildPreemptibleTestActorRef(attempt = 1, preemptible = 1, previousTransientRetriesCount = 0)
    val batchBackend = actorRef.underlyingActor
    val runId = generateStandardAsyncJob
    val handle = new GcpBatchPendingExecutionHandle(null, runId, None, None)

    val failedStatus = RunStatus.Failed(
      GcpBatchExitCode.VMPreemption,
      Seq.empty // task has not started
    )
    val executionResult = batchBackend.handleExecutionResult(failedStatus, handle)
    val result = Await.result(executionResult, timeout)
    result.isInstanceOf[FailedRetryableExecutionHandle] shouldBe true
    val failedHandle = result.asInstanceOf[FailedRetryableExecutionHandle]
    failedHandle.returnCode shouldBe None
    failedHandle.kvPairsToSave.map { pairs =>
      pairs.exists(p => p.key.key == GcpBatchBackendLifecycleActorFactory.preemptionCountKey && p.value == "0")
      pairs.exists(p => p.key.key == GcpBatchBackendLifecycleActorFactory.unexpectedRetryCountKey && p.value == "0")
      pairs.exists(p => p.key.key == GcpBatchBackendLifecycleActorFactory.transientRetryCountKey && p.value == "1")
    }
  }

  it should "choose preemptible over transient retry when task started before failing" in {
    val actorRef = buildPreemptibleTestActorRef(attempt = 1, preemptible = 1, previousTransientRetriesCount = 0)
    val batchBackend = actorRef.underlyingActor
    val runId = generateStandardAsyncJob
    val handle = new GcpBatchPendingExecutionHandle(null, runId, None, None)

    val failedStatus = RunStatus.Failed(
      GcpBatchExitCode.VMPreemption,
      List(ExecutionEvent(CallMetadataKeys.VmStartTime, OffsetDateTime.now()))
    )
    val executionResult = batchBackend.handleExecutionResult(failedStatus, handle)
    val result = Await.result(executionResult, timeout)
    result.isInstanceOf[FailedRetryableExecutionHandle] shouldBe true
    val failedHandle = result.asInstanceOf[FailedRetryableExecutionHandle]
    failedHandle.returnCode shouldBe None
    failedHandle.kvPairsToSave.map { pairs =>
      pairs.exists(p => p.key.key == GcpBatchBackendLifecycleActorFactory.preemptionCountKey && p.value == "2")
      pairs.exists(p => p.key.key == GcpBatchBackendLifecycleActorFactory.unexpectedRetryCountKey && p.value == "0")
      pairs.exists(p => p.key.key == GcpBatchBackendLifecycleActorFactory.transientRetryCountKey && p.value == "0")
    }
  }

  it should "handle Failure Status for various errors" in {

    val actorRef = buildPreemptibleTestActorRef(1, 1, 0)
    val batchBackend = actorRef.underlyingActor
    val runId = generateStandardAsyncJob
    val handle = new GcpBatchPendingExecutionHandle(null, runId, None, None)

    def checkFailedResult(errorCode: GcpBatchExitCode, events: List[ExecutionEvent] = List.empty): ExecutionHandle = {
      val failed = RunStatus.Failed(
        errorCode,
        events
      )
      Await.result(batchBackend.handleExecutionResult(failed, handle), timeout)
    }

    // Should retry
    checkFailedResult(GcpBatchExitCode.VMPreemption)
      .isInstanceOf[FailedRetryableExecutionHandle] shouldBe true
    checkFailedResult(GcpBatchExitCode.VMRecreatedDuringExecution) // no VM start time - task has not started
      .isInstanceOf[FailedRetryableExecutionHandle] shouldBe true

    // Should not retry
    checkFailedResult(GcpBatchExitCode.Success)
      .isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    checkFailedResult(GcpBatchExitCode.TaskRunsOverMaximumRuntime)
      .isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    checkFailedResult(
      GcpBatchExitCode.VMRecreatedDuringExecution,
      List(ExecutionEvent("Job state is set from SCHEDULED to RUNNING for job f00b4r", OffsetDateTime.now()))
    ).isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    actorRef.stop()
  }

  it should "map GCS paths and *only* GCS paths to local" in {
    val wdlString =
      s"""|workflow wf {
          |  call t
          |}
          |
          |task t {
          |  String abc
          |  File lf
          |  File gcsf
          |  command {}
          |  runtime { docker: "ubuntu" }
          |}
          |""".stripMargin

    val stringKey = "wf.t.abc"
    val stringVal = WomString("abc")
    val localFileKey = "wf.t.lf"
    val localFileVal = WomSingleFile("/blah/abc")
    val gcsFileKey = "wf.t.gcsf"
    val gcsFileVal = WomSingleFile("gs://blah/abc")

    val inputs: Map[String, WomValue] = Map(
      stringKey -> stringVal,
      localFileKey -> localFileVal,
      gcsFileKey -> gcsFileVal
    )

    val wdlNamespace = WdlNamespaceWithWorkflow
      .load(
        wdlString,
        Seq.empty[Draft2ImportResolver]
      )
      .get
    val womWorkflow = wdlNamespace.workflow
      .toWomWorkflowDefinition(isASubworkflow = false)
      .getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    wdlNamespace.toWomExecutable(Option(inputs.toJson.compactPrint), NoIoFunctionSet, strictValidation = true) match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap { case (port, v) =>
          v.select[WomValue] map { port -> _ }
        }

        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty,
          HogGroup("foo"),
          List.empty,
          None
        )

        val call: CommandCallNode = workflowDescriptor.callable.graph.nodes.collectFirst { case t: CommandCallNode =>
          t
        }.get
        val key = BackendJobDescriptorKey(call, None, 1)
        val runtimeAttributes = makeRuntimeAttributes(call)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor,
                                                 key,
                                                 runtimeAttributes,
                                                 fqnWdlMapToDeclarationMap(inputs),
                                                 NoDocker,
                                                 None,
                                                 Map.empty
        )

        val props = Props(new TestableGcpBatchJobExecutionActor(jobDescriptor, Promise(), gcpBatchConfiguration))
        val testActorRef = TestActorRef[TestableGcpBatchJobExecutionActor](
          props,
          s"TestableGcpBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
        )

        def gcsPathToLocal(womValue: WomValue): WomValue =
          WomFileMapper.mapWomFiles(testActorRef.underlyingActor.mapCommandLineWomFile)(womValue).get

        val mappedInputs = jobDescriptor.localInputs safeMapValues gcsPathToLocal

        mappedInputs(stringKey) match {
          case WomString(v) => assert(v.equalsIgnoreCase(stringVal.value))
          case _ => fail("test setup error")
        }

        mappedInputs(localFileKey) match {
          case wdlFile: WomSingleFile => assert(wdlFile.value.equalsIgnoreCase(localFileVal.value))
          case _ => fail("test setup error")
        }

        mappedInputs(gcsFileKey) match {
          case wdlFile: WomSingleFile => wdlFile.value shouldBe s"$MountPoint/blah/abc"
          case _ => fail("test setup error")
        }
      case Left(badtimes) => fail(badtimes.toList.mkString(", "))
    }
  }

  private val dockerAndDiskMapsWdlNamespace =
    WdlNamespaceWithWorkflow
      .load(
        SampleWdl.CurrentDirectoryMaps.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
        Seq.empty[Draft2ImportResolver]
      )
      .get

  private val dockerAndDiskArrayWdlNamespace =
    WdlNamespaceWithWorkflow
      .load(
        SampleWdl.CurrentDirectoryArray.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
        Seq.empty[Draft2ImportResolver]
      )
      .get

  private val dockerAndDiskFilesWdlNamespace =
    WdlNamespaceWithWorkflow
      .load(
        SampleWdl.CurrentDirectoryFiles.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
        Seq.empty[Draft2ImportResolver]
      )
      .get

  it should "generate correct BatchFileInputs from a WdlMap" in {
    val inputs: Map[String, WomValue] = Map(
      "stringToFileMap" -> WomMap(
        WomMapType(WomStringType, WomSingleFileType),
        Map(
          WomString("stringTofile1") -> WomSingleFile("gs://path/to/stringTofile1"),
          WomString("stringTofile2") -> WomSingleFile("gs://path/to/stringTofile2")
        )
      ),
      "fileToStringMap" -> WomMap(
        WomMapType(WomSingleFileType, WomStringType),
        Map(
          WomSingleFile("gs://path/to/fileToString1") -> WomString("fileToString1"),
          WomSingleFile("gs://path/to/fileToString2") -> WomString("fileToString2")
        )
      ),
      "fileToFileMap" -> WomMap(
        WomMapType(WomSingleFileType, WomSingleFileType),
        Map(
          WomSingleFile("gs://path/to/fileToFile1Key") -> WomSingleFile("gs://path/to/fileToFile1Value"),
          WomSingleFile("gs://path/to/fileToFile2Key") -> WomSingleFile("gs://path/to/fileToFile2Value")
        )
      ),
      "stringToString" -> WomMap(
        WomMapType(WomStringType, WomStringType),
        Map(
          WomString("stringToString1") -> WomString("path/to/stringToString1"),
          WomString("stringToString2") -> WomString("path/to/stringToString2")
        )
      )
    )

    val workflowInputs = inputs map { case (key, value) =>
      (s"wf_whereami.whereami.$key", value)
    }

    val womWorkflow =
      dockerAndDiskMapsWdlNamespace.workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    val womExecutableChecked =
      dockerAndDiskMapsWdlNamespace
        .toWomExecutable(Option(workflowInputs.toJson.compactPrint), NoIoFunctionSet, strictValidation = true)
    womExecutableChecked match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap { case (port, v) =>
          v.select[WomValue] map {
            port -> _
          }
        }
        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty,
          HogGroup("foo"),
          List.empty,
          None
        )

        val job: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
        val runtimeAttributes = makeRuntimeAttributes(job)
        val key = BackendJobDescriptorKey(job, None, 1)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor,
                                                 key,
                                                 runtimeAttributes,
                                                 fqnWdlMapToDeclarationMap(inputs),
                                                 NoDocker,
                                                 None,
                                                 Map.empty
        )

        val props = Props(new TestableGcpBatchJobExecutionActor(jobDescriptor, Promise(), gcpBatchConfiguration))
        val testActorRef = TestActorRef[TestableGcpBatchJobExecutionActor](
          props,
          s"TestableBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
        )

        val batchInputs = testActorRef.underlyingActor.generateInputs()
        batchInputs should have size 8
        batchInputs should contain(
          GcpBatchFileInput(
            name = "stringToFileMap",
            cloudPath = gcsPath("gs://path/to/stringTofile1"),
            relativeHostPath = DefaultPathBuilder.get("path/to/stringTofile1"),
            mount = workingDisk
          )
        )
        batchInputs should contain(
          GcpBatchFileInput(
            name = "stringToFileMap",
            cloudPath = gcsPath("gs://path/to/stringTofile2"),
            relativeHostPath = DefaultPathBuilder.get("path/to/stringTofile2"),
            mount = workingDisk
          )
        )
        batchInputs should contain(
          GcpBatchFileInput(
            name = "fileToStringMap",
            cloudPath = gcsPath("gs://path/to/fileToString1"),
            relativeHostPath = DefaultPathBuilder.get("path/to/fileToString1"),
            mount = workingDisk
          )
        )
        batchInputs should contain(
          GcpBatchFileInput(
            name = "fileToStringMap",
            cloudPath = gcsPath("gs://path/to/fileToString2"),
            relativeHostPath = DefaultPathBuilder.get("path/to/fileToString2"),
            mount = workingDisk
          )
        )
        batchInputs should contain(
          GcpBatchFileInput(
            name = "fileToFileMap",
            cloudPath = gcsPath("gs://path/to/fileToFile1Key"),
            relativeHostPath = DefaultPathBuilder.get("path/to/fileToFile1Key"),
            mount = workingDisk
          )
        )
        batchInputs should contain(
          GcpBatchFileInput(
            name = "fileToFileMap",
            cloudPath = gcsPath("gs://path/to/fileToFile1Value"),
            relativeHostPath = DefaultPathBuilder.get("path/to/fileToFile1Value"),
            mount = workingDisk
          )
        )
        batchInputs should contain(
          GcpBatchFileInput(
            name = "fileToFileMap",
            cloudPath = gcsPath("gs://path/to/fileToFile2Key"),
            relativeHostPath = DefaultPathBuilder.get("path/to/fileToFile2Key"),
            mount = workingDisk
          )
        )
        batchInputs should contain(
          GcpBatchFileInput(
            name = "fileToFileMap",
            cloudPath = gcsPath("gs://path/to/fileToFile2Value"),
            relativeHostPath = DefaultPathBuilder.get("path/to/fileToFile2Value"),
            mount = workingDisk
          )
        )

      case Left(badness) => fail(badness.toList.mkString(", "))
    }
  }

  private def makeBatchActorRef(sampleWdl: SampleWdl,
                                workflowInputs: Map[FullyQualifiedName, WomValue],
                                callName: LocallyQualifiedName,
                                callInputs: Map[LocallyQualifiedName, WomValue],
                                functions: BatchExpressionFunctions = TestableGcpBatchExpressionFunctions
  ): TestActorRef[TestableGcpBatchJobExecutionActor] = {
    val wdlNamespaceWithWorkflow =
      WdlNamespaceWithWorkflow
        .load(
          sampleWdl.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
          Seq.empty[Draft2ImportResolver]
        )
        .get
    val womWorkflow =
      wdlNamespaceWithWorkflow.workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    val womExecutableChecked =
      wdlNamespaceWithWorkflow
        .toWomExecutable(
          Option(workflowInputs.toJson.compactPrint),
          NoIoFunctionSet,
          strictValidation = true
        )
    womExecutableChecked match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap { case (port, v) =>
          v.select[WomValue] map {
            port -> _
          }
        }
        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty,
          HogGroup("foo"),
          List.empty,
          None
        )

        val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.find(_.localName == callName).get
        val key = BackendJobDescriptorKey(call, None, 1)
        val runtimeAttributes = makeRuntimeAttributes(call)
        val jobDescriptor =
          BackendJobDescriptor(
            workflowDescriptor = workflowDescriptor,
            key = key,
            runtimeAttributes = runtimeAttributes,
            evaluatedTaskInputs = fqnWdlMapToDeclarationMap(callInputs),
            maybeCallCachingEligible = NoDocker,
            dockerSize = None,
            prefetchedKvStoreEntries = Map.empty
          )

        val props = Props(
          new TestableGcpBatchJobExecutionActor(jobDescriptor, Promise(), gcpBatchConfiguration, functions)
        )
        TestActorRef[TestableGcpBatchJobExecutionActor](
          props,
          s"TestableBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
        )
      case Left(badness) => fail(badness.toList.mkString(", "))
    }
  }

  it should "generate correct BatchOutputs" in {
    val womFile = WomSingleFile("gs://blah/b/c.txt")
    val workflowInputs = Map("file_passing.f" -> womFile)
    val callInputs = Map(
      "in" -> womFile, // how does one programmatically map the wf inputs to the call inputs?
      "out_name" -> WomString("out") // is it expected that this isn't using the default?
    )
    val batchBackend = makeBatchActorRef(SampleWdl.FilePassingWorkflow, workflowInputs, "a", callInputs).underlyingActor
    val jobDescriptor = batchBackend.jobDescriptor
    val workflowId = batchBackend.workflowId
    val batchInputs = batchBackend.generateInputs()
    batchInputs should have size 1
    batchInputs should contain(
      GcpBatchFileInput(
        name = "in",
        cloudPath = gcsPath("gs://blah/b/c.txt"),
        relativeHostPath = DefaultPathBuilder.get("blah/b/c.txt"),
        mount = workingDisk
      )
    )
    val batchOutputs = batchBackend.generateOutputs(jobDescriptor)
    batchOutputs should have size 1
    batchOutputs should contain(
      GcpBatchFileOutput(
        "out",
        gcsPath(s"gs://my-cromwell-workflows-bucket/file_passing/$workflowId/call-a/out"),
        DefaultPathBuilder.get("out"),
        workingDisk,
        optional = false,
        secondary = false
      )
    )
  }

  it should "generate correct BatchInputs when a command line contains a write_lines call in it" in {
    val inputs = Map(
      "strs" -> WomArray(WomArrayType(WomStringType), Seq("A", "B", "C").map(WomString))
    )

    class TestBatchApiExpressionFunctions extends BatchExpressionFunctions(TestableStandardExpressionFunctionsParams) {
      override def writeFile(path: String, content: String): Future[WomSingleFile] =
        Future.fromTry(Success(WomSingleFile(s"gs://some/path/file.txt")))
    }

    val functions = new TestBatchApiExpressionFunctions
    val batchBackend = makeBatchActorRef(SampleWdl.ArrayIO, Map.empty, "serialize", inputs, functions).underlyingActor
    val jobDescriptor = batchBackend.jobDescriptor
    val batchInputs = batchBackend.generateInputs()
    batchInputs should have size 1
    batchInputs should contain(
      GcpBatchFileInput(
        name = "c35ad8d3",
        cloudPath = gcsPath("gs://some/path/file.txt"),
        relativeHostPath = DefaultPathBuilder.get("some/path/file.txt"),
        mount = workingDisk
      )
    )
    val batchOutputs = batchBackend.generateOutputs(jobDescriptor)
    batchOutputs should have size 0
  }

  it should "generate correct BatchFileInputs from a WdlArray" in {
    val inputs: Map[String, WomValue] = Map(
      "fileArray" ->
        WomArray(WomArrayType(WomSingleFileType),
                 Seq(WomSingleFile("gs://path/to/file1"), WomSingleFile("gs://path/to/file2"))
        )
    )

    val workflowInputs = inputs map { case (key, value) =>
      (s"wf_whereami.whereami.$key", value)
    }

    val womWorkflow =
      dockerAndDiskArrayWdlNamespace.workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    val womExecutableChecked =
      dockerAndDiskArrayWdlNamespace
        .toWomExecutable(Option(workflowInputs.toJson.compactPrint), NoIoFunctionSet, strictValidation = true)
    womExecutableChecked match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap { case (port, v) =>
          v.select[WomValue] map {
            port -> _
          }
        }
        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty,
          HogGroup("foo"),
          List.empty,
          None
        )

        val job: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
        val runtimeAttributes = makeRuntimeAttributes(job)
        val key = BackendJobDescriptorKey(job, None, 1)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor,
                                                 key,
                                                 runtimeAttributes,
                                                 fqnWdlMapToDeclarationMap(inputs),
                                                 NoDocker,
                                                 None,
                                                 Map.empty
        )

        val props = Props(new TestableGcpBatchJobExecutionActor(jobDescriptor, Promise(), gcpBatchConfiguration))
        val testActorRef = TestActorRef[TestableGcpBatchJobExecutionActor](
          props,
          s"TestableBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
        )

        val batchInputs = testActorRef.underlyingActor.generateInputs()
        batchInputs should have size 2
        batchInputs should contain(
          GcpBatchFileInput(
            name = "fileArray",
            cloudPath = gcsPath("gs://path/to/file1"),
            relativeHostPath = DefaultPathBuilder.get("path/to/file1"),
            mount = workingDisk
          )
        )
        batchInputs should contain(
          GcpBatchFileInput(
            name = "fileArray",
            cloudPath = gcsPath("gs://path/to/file2"),
            relativeHostPath = DefaultPathBuilder.get("path/to/file2"),
            mount = workingDisk
          )
        )
      case Left(badness) => fail(badness.toList.mkString(", "))
    }
  }

  it should "generate correct BatchFileInputs from a WdlFile" in {
    val inputs: Map[String, WomValue] = Map(
      "file1" -> WomSingleFile("gs://path/to/file1"),
      "file2" -> WomSingleFile("gs://path/to/file2")
    )

    val workflowInputs = inputs map { case (key, value) =>
      (s"wf_whereami.whereami.$key", value)
    }

    val womWorkflow =
      dockerAndDiskFilesWdlNamespace.workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    val womExecutableChecked =
      dockerAndDiskFilesWdlNamespace
        .toWomExecutable(Option(workflowInputs.toJson.compactPrint), NoIoFunctionSet, strictValidation = true)
    womExecutableChecked match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap { case (port, v) =>
          v.select[WomValue] map {
            port -> _
          }
        }
        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty,
          HogGroup("foo"),
          List.empty,
          None
        )

        val job: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
        val runtimeAttributes = makeRuntimeAttributes(job)
        val key = BackendJobDescriptorKey(job, None, 1)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor,
                                                 key,
                                                 runtimeAttributes,
                                                 fqnWdlMapToDeclarationMap(inputs),
                                                 NoDocker,
                                                 None,
                                                 Map.empty
        )

        val props = Props(new TestableGcpBatchJobExecutionActor(jobDescriptor, Promise(), gcpBatchConfiguration))
        val testActorRef = TestActorRef[TestableGcpBatchJobExecutionActor](
          props,
          s"TestableBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
        )

        val batchInputs = testActorRef.underlyingActor.generateInputs()

        batchInputs should have size 2
        batchInputs should contain(
          GcpBatchFileInput(
            name = "file1",
            cloudPath = gcsPath("gs://path/to/file1"),
            relativeHostPath = DefaultPathBuilder.get("path/to/file1"),
            mount = workingDisk
          )
        )
        batchInputs should contain(
          GcpBatchFileInput(
            name = "file2",
            cloudPath = gcsPath("gs://path/to/file2"),
            relativeHostPath = DefaultPathBuilder.get("path/to/file2"),
            mount = workingDisk
          )
        )

      case Left(badness) => fail(badness.toList.mkString(", "))
    }
  }

  it should "convert local Paths back to corresponding GCS paths in BatchOutputs" in {

    val batchOutputs = Set(
      GcpBatchFileOutput(
        s"$MountPoint/centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file1",
        gcsPath("gs://centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file1"),
        DefaultPathBuilder.get(s"$MountPoint/centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file1"),
        workingDisk,
        optional = false,
        secondary = false
      ),
      GcpBatchFileOutput(
        s"$MountPoint/centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file2",
        gcsPath("gs://centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file2"),
        DefaultPathBuilder.get(s"$MountPoint/centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file2"),
        workingDisk,
        optional = false,
        secondary = false
      ),
      GcpBatchFileOutput(
        s"$MountPoint/centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file3",
        gcsPath("gs://centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file3"),
        DefaultPathBuilder.get(s"$MountPoint/centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file3"),
        workingDisk,
        optional = false,
        secondary = false
      ),
      GcpBatchFileOutput(
        s"$MountPoint/centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file4",
        gcsPath("gs://centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file4"),
        DefaultPathBuilder.get(s"$MountPoint/centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file4"),
        workingDisk,
        optional = false,
        secondary = false
      ),
      GcpBatchFileOutput(
        s"$MountPoint/centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file5",
        gcsPath("gs://centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file5"),
        DefaultPathBuilder.get(s"$MountPoint/centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file5"),
        workingDisk,
        optional = false,
        secondary = false
      )
    )
    val outputValues = Seq(
      WomSingleFile(s"$MountPoint/centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file1"),
      WomArray(
        WomArrayType(WomSingleFileType),
        Seq(
          WomSingleFile(s"$MountPoint/centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file2"),
          WomSingleFile(s"$MountPoint/centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file3")
        )
      ),
      WomMap(
        WomMapType(WomSingleFileType, WomSingleFileType),
        Map(
          WomSingleFile(
            s"$MountPoint/centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file4"
          ) -> WomSingleFile(s"$MountPoint/centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file5")
        )
      )
    )

    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      WdlNamespaceWithWorkflow
        .load(SampleWdl.EmptyString.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
              Seq.empty[Draft2ImportResolver]
        )
        .get
        .workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      NoOptions,
      Labels.empty,
      HogGroup("foo"),
      List.empty,
      None
    )

    val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
    val key = BackendJobDescriptorKey(call, None, 1)
    val runtimeAttributes = makeRuntimeAttributes(call)
    val jobDescriptor =
      BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, None, Map.empty)

    val props = Props(new TestableGcpBatchJobExecutionActor(jobDescriptor, Promise(), gcpBatchConfiguration))
    val testActorRef = TestActorRef[TestableGcpBatchJobExecutionActor](
      props,
      s"TestableBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
    )

    def wdlValueToGcsPath(batchOutputs: Set[GcpBatchFileOutput])(womValue: WomValue): WomValue =
      WomFileMapper.mapWomFiles(testActorRef.underlyingActor.womFileToGcsPath(batchOutputs.toSet))(womValue).get

    val result = outputValues map wdlValueToGcsPath(batchOutputs)
    result should have size 3
    result should contain(WomSingleFile("gs://centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file1"))
    result should contain(
      WomArray(
        WomArrayType(WomSingleFileType),
        Seq(
          WomSingleFile("gs://centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file2"),
          WomSingleFile("gs://centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file3")
        )
      )
    )
    result should contain(
      WomMap(
        WomMapType(WomSingleFileType, WomSingleFileType),
        Map(
          WomSingleFile("gs://centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file4") -> WomSingleFile(
            "gs://centaur-ci-public/GcpBatchAsyncBackendJobExecutionActorSpec/file5"
          )
        )
      )
    )
  }

  it should "create a GcpBatchFileInput for the monitoring script, when specified" in {

    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      WdlNamespaceWithWorkflow
        .load(SampleWdl.EmptyString.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
              Seq.empty[Draft2ImportResolver]
        )
        .get
        .workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      WorkflowOptions.fromJsonString("""{"monitoring_script": "gs://path/to/script"}""").get,
      Labels.empty,
      HogGroup("foo"),
      List.empty,
      None
    )

    val job: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
    val runtimeAttributes = makeRuntimeAttributes(job)
    val key = BackendJobDescriptorKey(job, None, 1)
    val jobDescriptor =
      BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, None, Map.empty)

    val props = Props(new TestableGcpBatchJobExecutionActor(jobDescriptor, Promise(), gcpBatchConfiguration))
    val testActorRef = TestActorRef[TestableGcpBatchJobExecutionActor](
      props,
      s"TestableGcpBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
    )

    testActorRef.underlyingActor.monitoringScript shouldBe
      Option(
        GcpBatchFileInput("monitoring-in",
                          gcsPath("gs://path/to/script"),
                          DefaultPathBuilder.get("monitoring.sh"),
                          workingDisk
        )
      )
  }

  it should "not create a GcpBatchFileInput for the monitoring script, when not specified" in {

    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      WdlNamespaceWithWorkflow
        .load(SampleWdl.EmptyString.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
              Seq.empty[Draft2ImportResolver]
        )
        .get
        .workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      NoOptions,
      Labels.empty,
      HogGroup("foo"),
      List.empty,
      None
    )

    val job: CommandCallNode = workflowDescriptor.callable.graph.nodes.collectFirst { case t: CommandCallNode => t }.get
    val key = BackendJobDescriptorKey(job, None, 1)
    val runtimeAttributes = makeRuntimeAttributes(job)
    val jobDescriptor =
      BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, None, Map.empty)

    val props = Props(new TestableGcpBatchJobExecutionActor(jobDescriptor, Promise(), gcpBatchConfiguration))
    val testActorRef = TestActorRef[TestableGcpBatchJobExecutionActor](
      props,
      s"TestableGcpBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
    )

    testActorRef.underlyingActor.monitoringScript shouldBe None
  }

  it should "return GCP Batch log paths for non-scattered call" in {

    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId(UUID.fromString("e6236763-c518-41d0-9688-432549a8bf7c")),
      WdlNamespaceWithWorkflow
        .load(SampleWdl.HelloWorld.asWorkflowSources(""" runtime {docker: "ubuntu:latest"} """).workflowSource.get,
              Seq.empty[Draft2ImportResolver]
        )
        .get
        .workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      WorkflowOptions
        .fromJsonString(s""" {"${GcpBatchWorkflowPaths.GcsRootOptionKey}": "gs://path/to/gcs_root"} """)
        .get,
      Labels.empty,
      HogGroup("foo"),
      List.empty,
      None
    )

    val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.find(_.localName == "hello").get
    val key = BackendJobDescriptorKey(call, None, 1)
    val runtimeAttributes = makeRuntimeAttributes(call)
    val jobDescriptor =
      BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, None, Map.empty)

    val props = Props(new TestableGcpBatchJobExecutionActor(jobDescriptor, Promise(), gcpBatchConfiguration))
    val testActorRef = TestActorRef[TestableGcpBatchJobExecutionActor](
      props,
      s"TestableGcpBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
    )

    val batchBackend = testActorRef.underlyingActor

    batchBackend.gcpBatchCallPaths.stdout should be(a[GcsPath])
    batchBackend.gcpBatchCallPaths.stdout.pathAsString shouldBe
      "gs://path/to/gcs_root/wf_hello/e6236763-c518-41d0-9688-432549a8bf7c/call-hello/stdout"
    batchBackend.gcpBatchCallPaths.stderr should be(a[GcsPath])
    batchBackend.gcpBatchCallPaths.stderr.pathAsString shouldBe
      "gs://path/to/gcs_root/wf_hello/e6236763-c518-41d0-9688-432549a8bf7c/call-hello/stderr"
    batchBackend.gcpBatchCallPaths.batchLogPath should be(a[GcsPath])
    batchBackend.gcpBatchCallPaths.batchLogPath.pathAsString shouldBe
      "gs://path/to/gcs_root/wf_hello/e6236763-c518-41d0-9688-432549a8bf7c/call-hello/hello.log"
  }

  it should "return Batch log paths for scattered call" in {

    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId(UUID.fromString("e6236763-c518-41d0-9688-432549a8bf7d")),
      WdlNamespaceWithWorkflow
        .load(
          new SampleWdl.ScatterWdl().asWorkflowSources(""" runtime {docker: "ubuntu:latest"} """).workflowSource.get,
          Seq.empty[Draft2ImportResolver]
        )
        .get
        .workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      WorkflowOptions
        .fromJsonString(s""" {"${GcpBatchWorkflowPaths.GcsRootOptionKey}": "gs://path/to/gcs_root"} """)
        .get,
      Labels.empty,
      HogGroup("foo"),
      List.empty,
      None
    )

    val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.find(_.localName == "B").get
    val key = BackendJobDescriptorKey(call, Option(2), 1)
    val runtimeAttributes = makeRuntimeAttributes(call)
    val jobDescriptor =
      BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, None, Map.empty)

    val props = Props(new TestableGcpBatchJobExecutionActor(jobDescriptor, Promise(), gcpBatchConfiguration))
    val testActorRef = TestActorRef[TestableGcpBatchJobExecutionActor](
      props,
      s"TestableGcpBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
    )

    val batchBackend = testActorRef.underlyingActor

    batchBackend.gcpBatchCallPaths.stdout should be(a[GcsPath])
    batchBackend.gcpBatchCallPaths.stdout.pathAsString shouldBe
      "gs://path/to/gcs_root/w/e6236763-c518-41d0-9688-432549a8bf7d/call-B/shard-2/stdout"
    batchBackend.gcpBatchCallPaths.stderr should be(a[GcsPath])
    batchBackend.gcpBatchCallPaths.stderr.pathAsString shouldBe
      "gs://path/to/gcs_root/w/e6236763-c518-41d0-9688-432549a8bf7d/call-B/shard-2/stderr"
    batchBackend.gcpBatchCallPaths.batchLogPath should be(a[GcsPath])
    batchBackend.gcpBatchCallPaths.batchLogPath.pathAsString shouldBe
      "gs://path/to/gcs_root/w/e6236763-c518-41d0-9688-432549a8bf7d/call-B/shard-2/B-2.log"
  }

  it should "return preemptible = true only in the correct cases" in {

    def attempt(max: Int, attempt: Int): GcpBatchAsyncBackendJobExecutionActor =
      buildPreemptibleTestActorRef(attempt, max, 0).underlyingActor
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

  it should "return the project from the workflow options in the start metadata" in {

    val googleProject = "baa-ram-ewe"
    val batchGcsRoot = "gs://anorexic/duck"
    val workflowId = WorkflowId.randomId()
    val workflowDescriptor = BackendWorkflowDescriptor(
      workflowId,
      WdlNamespaceWithWorkflow
        .load(
          SampleWdl.EmptyString.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
          Seq.empty[Draft2ImportResolver]
        )
        .get
        .workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      WorkflowOptions
        .fromJsonString(
          s"""|{
              |  "google_project": "$googleProject",
              |  "${GcpBatchWorkflowPaths.GcsRootOptionKey}": "$batchGcsRoot"
              |}
              |""".stripMargin
        )
        .get,
      Labels.empty,
      HogGroup("foo"),
      List.empty,
      None
    )

    val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.find(_.localName == "goodbye").get
    val key = BackendJobDescriptorKey(call, None, 1)
    val runtimeAttributes = makeRuntimeAttributes(call)
    val jobDescriptor =
      BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, None, Map.empty)

    val props = Props(new TestableGcpBatchJobExecutionActor(jobDescriptor, Promise(), gcpBatchConfiguration))
    val testActorRef = TestActorRef[TestableGcpBatchJobExecutionActor](
      props,
      s"TestableGcpBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
    )

    val batchBackend = testActorRef.underlyingActor

    val actual = batchBackend.startMetadataKeyValues.safeMapValues(_.toString)
    actual should be(
      Map(
        "callRoot" -> s"$batchGcsRoot/wf_hello/$workflowId/call-goodbye",
        "gcpBatch:executionBucket" -> batchGcsRoot,
        "gcpBatch:googleProject" -> googleProject,
        "labels:cromwell-workflow-id" -> s"cromwell-$workflowId",
        "labels:wdl-task-name" -> "goodbye",
        "preemptible" -> "false",
        "runtimeAttributes:bootDiskSizeGb" -> "30",
        "runtimeAttributes:continueOnReturnCode" -> "0",
        "runtimeAttributes:cpu" -> "1",
        "runtimeAttributes:disks" -> "local-disk 200 SSD",
        "runtimeAttributes:docker" -> "ubuntu:latest",
        "runtimeAttributes:failOnStderr" -> "false",
        "runtimeAttributes:memory" -> "2 GB",
        "runtimeAttributes:noAddress" -> "false",
        "runtimeAttributes:preemptible" -> "0",
        "runtimeAttributes:zones" -> "us-central1-b,us-central1-a",
        "runtimeAttributes:maxRetries" -> "0",
        "stderr" -> s"$batchGcsRoot/wf_hello/$workflowId/call-goodbye/stderr",
        "stdout" -> s"$batchGcsRoot/wf_hello/$workflowId/call-goodbye/stdout"
      )
    )
  }

  private def buildJobDescriptor(): BackendJobDescriptor = {
    val attempt = 1
    val preemptible = 1
    val wdlNamespace = WdlNamespaceWithWorkflow
      .load(YoSup.replace("[PREEMPTIBLE]", s"preemptible: $preemptible"), Seq.empty[Draft2ImportResolver])
      .get
    val womDefinition = wdlNamespace.workflow
      .toWomWorkflowDefinition(isASubworkflow = false)
      .getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))

    wdlNamespace.toWomExecutable(Option(Inputs.toJson.compactPrint), NoIoFunctionSet, strictValidation = true) match {
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
          Labels.empty,
          HogGroup("foo"),
          List.empty,
          None
        )
        val job = workflowDescriptor.callable.taskCallNodes.head
        val key = BackendJobDescriptorKey(job, None, attempt)
        val runtimeAttributes = makeRuntimeAttributes(job)

        BackendJobDescriptor(workflowDescriptor,
                             key,
                             runtimeAttributes,
                             fqnWdlMapToDeclarationMap(Inputs),
                             NoDocker,
                             None,
                             Map()
        )
      case Left(badtimes) => fail(badtimes.toList.mkString(", "))
    }
  }

  def buildTestActorRef(
    jobDescriptor: BackendJobDescriptor,
    serviceRegistryActor: Option[TestProbe] = None
  ): TestActorRef[TestableGcpBatchJobExecutionActor] = {
    val props = Props(
      new TestableGcpBatchJobExecutionActor(
        jobDescriptor,
        Promise(),
        gcpBatchConfiguration,
        TestableGcpBatchExpressionFunctions,
        emptyActor,
        failIoActor,
        serviceRegistryActor.map(actor => actor.ref).getOrElse(kvService)
      )
    )
    TestActorRef(props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")
  }

  it should "emit expected timing metadata as task executes" in {
    val expectedJobStart = OffsetDateTime.now().minus(3, ChronoUnit.HOURS)
    val expectedVmStart = OffsetDateTime.now().minus(2, ChronoUnit.HOURS)
    val expectedVmEnd = OffsetDateTime.now().minus(1, ChronoUnit.HOURS)

    val pollResult0 = RunStatus.Initializing(Seq.empty)
    val pollResult1 = RunStatus.Running(Seq(ExecutionEvent("fakeEvent", expectedJobStart)))
    val pollResult2 = RunStatus.Running(Seq(ExecutionEvent(CallMetadataKeys.VmStartTime, expectedVmStart)))
    val pollResult3 = RunStatus.Running(Seq(ExecutionEvent(CallMetadataKeys.VmEndTime, expectedVmEnd)))
    val terminalPollResult =
      RunStatus.Success(Seq(ExecutionEvent("fakeEvent", OffsetDateTime.now().truncatedTo(ChronoUnit.MILLIS))))

    val serviceRegistryProbe = TestProbe()

    val jobDescriptor = buildJobDescriptor()
    val job = StandardAsyncJob(UUID.randomUUID().toString)
    val run = Run(job)
    val handle = new GcpBatchPendingExecutionHandle(jobDescriptor, run.job, Option(run), None)
    val testActorRef = buildTestActorRef(jobDescriptor, Option(serviceRegistryProbe))

    testActorRef.underlyingActor.handlePollSuccess(handle, pollResult0)
    serviceRegistryProbe.fishForMessage(max = 5.seconds.dilated, hint = "") {
      case _: PutMetadataAction => true
      case _ => false
    }
    testActorRef.underlyingActor.handlePollSuccess(handle, pollResult1)
    serviceRegistryProbe.fishForMessage(max = 5.seconds.dilated, hint = "") {
      case _: PutMetadataAction => true
      case _ => false
    }
    testActorRef.underlyingActor.handlePollSuccess(handle, pollResult2)
    serviceRegistryProbe.fishForMessage(max = 5.seconds.dilated, hint = "") {
      case action: PutMetadataAction =>
        action.events.exists(event => event.key.key.equals(CallMetadataKeys.VmStartTime))
      case _ => false
    }
    testActorRef.underlyingActor.handlePollSuccess(handle, pollResult3)
    serviceRegistryProbe.fishForMessage(max = 5.seconds.dilated, hint = "") {
      case action: PutMetadataAction => action.events.exists(event => event.key.key.equals(CallMetadataKeys.VmEndTime))
      case _ => false
    }
    testActorRef.underlyingActor.handlePollSuccess(handle, terminalPollResult)
    serviceRegistryProbe.fishForMessage(max = 5.seconds.dilated, hint = "") {
      case _: PutMetadataAction => true
      case _ => false
    }
  }

  it should "send bard metrics message on task success" in {
    val expectedJobStart = OffsetDateTime.now().minus(3, ChronoUnit.HOURS)
    val expectedVmStart = OffsetDateTime.now().minus(2, ChronoUnit.HOURS)
    val expectedVmEnd = OffsetDateTime.now().minus(1, ChronoUnit.HOURS)

    val pollResult0 = RunStatus.Initializing(Seq.empty)
    val pollResult1 = RunStatus.Running(Seq(ExecutionEvent("fakeEvent", expectedJobStart)))
    val pollResult2 = RunStatus.Running(Seq(ExecutionEvent(CallMetadataKeys.VmStartTime, expectedVmStart)))
    val pollResult3 = RunStatus.Running(Seq(ExecutionEvent(CallMetadataKeys.VmEndTime, expectedVmEnd)))
    val terminalPollResult =
      RunStatus.Success(Seq(ExecutionEvent("fakeEvent", OffsetDateTime.now().truncatedTo(ChronoUnit.MILLIS))))

    val serviceRegistryProbe = TestProbe()

    val jobDescriptor = buildJobDescriptor()
    val job = StandardAsyncJob(jobDescriptor.workflowDescriptor.id.id.toString)
    val run = Run(job)
    val handle = new GcpBatchPendingExecutionHandle(jobDescriptor, run.job, Option(run), None)
    val testActorRef = buildTestActorRef(jobDescriptor, Option(serviceRegistryProbe))

    testActorRef.underlyingActor.handlePollSuccess(handle, pollResult0)
    testActorRef.underlyingActor.handlePollSuccess(handle, pollResult1)
    testActorRef.underlyingActor.handlePollSuccess(handle, pollResult2)
    testActorRef.underlyingActor.handlePollSuccess(handle, pollResult3)
    testActorRef.underlyingActor.handlePollSuccess(handle, terminalPollResult)

    val bardMessage = serviceRegistryProbe.fishForMessage(5.seconds) {
      case _: BardEventRequest => true
      case _ => false
    }
    val taskSummary = bardMessage.asInstanceOf[BardEventRequest].event.asInstanceOf[TaskSummaryEvent]
    taskSummary.workflowId should be(jobDescriptor.workflowDescriptor.id.id)
    taskSummary.parentWorkflowId should be(None)
    taskSummary.rootWorkflowId should be(jobDescriptor.workflowDescriptor.id.id)
    taskSummary.jobTag should be(jobDescriptor.key.tag)
    taskSummary.jobFullyQualifiedName should be(jobDescriptor.key.call.fullyQualifiedName)
    taskSummary.jobIndex should be(None)
    taskSummary.jobAttempt should be(jobDescriptor.key.attempt)
    taskSummary.terminalState shouldBe a[String]
    taskSummary.platform should be(Some("gcp"))
    taskSummary.dockerImage should be(Some("ubuntu:latest"))
    taskSummary.cpuCount should be(1)
    taskSummary.memoryBytes should be(2.147483648e9)
    taskSummary.startTime should not be empty
    taskSummary.cpuStartTime should be(Option(expectedVmStart.toString))
    taskSummary.endTime should not be empty
    taskSummary.jobSeconds should be(7200)
    taskSummary.cpuSeconds should be(Option(3600))
  }

  it should "send bard metrics message on task failure" in {
    val expectedJobStart = OffsetDateTime.now().minus(3, ChronoUnit.HOURS)
    val expectedVmStart = OffsetDateTime.now().minus(2, ChronoUnit.HOURS)

    val pollResult0 = RunStatus.Initializing(Seq.empty)
    val pollResult1 = RunStatus.Running(Seq(ExecutionEvent("fakeEvent", expectedJobStart)))
    val pollResult2 = RunStatus.Running(Seq(ExecutionEvent(CallMetadataKeys.VmStartTime, expectedVmStart)))
    val abortStatus = RunStatus.Aborted()

    val serviceRegistryProbe = TestProbe()

    val jobDescriptor = buildJobDescriptor()
    val job = StandardAsyncJob(jobDescriptor.workflowDescriptor.id.id.toString)
    val run = Run(job)
    val handle = new GcpBatchPendingExecutionHandle(jobDescriptor, run.job, Option(run), None)
    val testActorRef = buildTestActorRef(jobDescriptor, Option(serviceRegistryProbe))

    testActorRef.underlyingActor.handlePollSuccess(handle, pollResult0)
    testActorRef.underlyingActor.handlePollSuccess(handle, pollResult1)
    testActorRef.underlyingActor.handlePollSuccess(handle, pollResult2)
    testActorRef.underlyingActor.handlePollSuccess(handle, abortStatus)

    val bardMessage = serviceRegistryProbe.fishForMessage(5.seconds) {
      case _: BardEventRequest => true
      case _ => false
    }

    val taskSummary = bardMessage.asInstanceOf[BardEventRequest].event.asInstanceOf[TaskSummaryEvent]
    taskSummary.workflowId should be(jobDescriptor.workflowDescriptor.id.id)
    taskSummary.parentWorkflowId should be(None)
    taskSummary.rootWorkflowId should be(jobDescriptor.workflowDescriptor.id.id)
    taskSummary.jobTag should be(jobDescriptor.key.tag)
    taskSummary.jobFullyQualifiedName should be(jobDescriptor.key.call.fullyQualifiedName)
    taskSummary.jobIndex should be(None)
    taskSummary.jobAttempt should be(jobDescriptor.key.attempt)
    taskSummary.terminalState shouldBe a[String]
    taskSummary.platform should be(Some("gcp"))
    taskSummary.dockerImage should be(Some("ubuntu:latest"))
    taskSummary.cpuCount should be(1)
    taskSummary.memoryBytes should be(2.147483648e9)
    taskSummary.startTime should not be empty
    taskSummary.cpuStartTime should be(Option(expectedVmStart.toString))
    taskSummary.endTime should not be empty
    taskSummary.jobSeconds should be > 0.toLong
    taskSummary.cpuSeconds.get should be > 0.toLong
  }

  private def makeRuntimeAttributes(job: CommandCallNode) = {
    val evaluatedAttributes =
      RuntimeAttributeDefinition.evaluateRuntimeAttributes(job.callable.runtimeAttributes, null, Map.empty)
    RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributesBuilder.definitions.toSet, NoOptions)(
      evaluatedAttributes.getOrElse(fail("Failed to evaluate runtime attributes"))
    )
  }

  private def generateStandardAsyncJob =
    StandardAsyncJob(
      JobName.newBuilder().setJob(UUID.randomUUID().toString).setProject("test").setLocation("local").build().toString
    )
}
