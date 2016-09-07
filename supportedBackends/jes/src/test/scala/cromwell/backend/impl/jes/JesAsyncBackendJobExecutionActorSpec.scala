package cromwell.backend.impl.jes

import java.nio.file.Paths
import java.util.UUID

import akka.actor.{ActorRef, Props}
import akka.event.LoggingAdapter
import akka.testkit.{ImplicitSender, TestActorRef, TestDuration}
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.async.AsyncBackendJobExecutionActor.{Execute, ExecutionMode}
import cromwell.backend.async.{AbortedExecutionHandle, ExecutionHandle, FailedNonRetryableExecutionHandle, FailedRetryableExecutionHandle}
import cromwell.backend.impl.jes.JesAsyncBackendJobExecutionActor.JesPendingExecutionHandle
import cromwell.backend.impl.jes.RunStatus.Failed
import cromwell.backend.impl.jes.io.{DiskType, JesWorkingDisk}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobDescriptorKey, BackendWorkflowDescriptor, PreemptedException, RuntimeAttributeDefinition}
import cromwell.core._
import cromwell.core.logging.LoggerWrapper
import cromwell.filesystems.gcs._
import cromwell.util.SampleWdl
import org.scalatest._
import cromwell.core.{WorkflowId, WorkflowOptions}
import cromwell.filesystems.gcs.GoogleAuthMode.GoogleAuthOptions
import org.scalatest.prop.Tables.Table
import org.slf4j.Logger
import org.specs2.mock.Mockito
import spray.json.{JsObject, JsValue}
import wdl4s.types.{WdlArrayType, WdlFileType, WdlMapType, WdlStringType}
import wdl4s.values.{WdlArray, WdlFile, WdlMap, WdlString, WdlValue}
import wdl4s.{Call, LocallyQualifiedName, NamespaceWithWorkflow}

import scala.concurrent.duration._
import scala.concurrent.{Await, ExecutionContext, Future, Promise}
import scala.util.{Success, Try}

class JesAsyncBackendJobExecutionActorSpec extends TestKitSuite("JesAsyncBackendJobExecutionActorSpec")
  with FlatSpecLike with Matchers with ImplicitSender with Mockito {

  import JesTestConfig._

  implicit val Timeout = 5.seconds.dilated

  val YoSup =
    """
      |task sup {
      |  String addressee
      |  command {
      |    echo "yo sup ${addressee}!"
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
      |workflow sup {
      |  call sup
      |}
    """.stripMargin

  val Inputs = Map("sup.sup.addressee" -> WdlString("dog"))

  val NoOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue]))

  val TestableCallContext = CallContext(MockGcsFileSystemBuilder.mockGcsFileSystem.getPath("gs://root"), "out", "err")

  val TestableJesExpressionFunctions = {
    new JesExpressionFunctions(List(MockGcsFileSystemBuilder.mockGcsFileSystem), TestableCallContext)
  }

  private def buildInitializationData(jobDescriptor: BackendJobDescriptor, configuration: JesConfiguration) = {
    def gcsFileSystem = {
      val authOptions = new GoogleAuthOptions {
        override def get(key: String): Try[String] = Try(throw new RuntimeException(s"key '$key' not found"))
      }

      val storage = jesConfiguration.jesAttributes.gcsFilesystemAuth.buildStorage(authOptions, jesConfiguration.googleConfig)
      GcsFileSystem(GcsFileSystemProvider(storage)(scala.concurrent.ExecutionContext.global))
    }

    val workflowPaths = JesWorkflowPaths(jobDescriptor.workflowDescriptor, configuration, gcsFileSystem)
    JesBackendInitializationData(workflowPaths, null)
  }

  class TestableJesJobExecutionActor(jobDescriptor: BackendJobDescriptor,
                                     promise: Promise[BackendJobExecutionResponse],
                                     jesConfiguration: JesConfiguration,
                                     functions: JesExpressionFunctions = TestableJesExpressionFunctions)
    extends JesAsyncBackendJobExecutionActor(jobDescriptor, promise, jesConfiguration, buildInitializationData(jobDescriptor, jesConfiguration), emptyActor) {

    override lazy val jobLogger = new LoggerWrapper {
      override def akkaLogger: Option[LoggingAdapter] = Option(log)

      override def tag: String = s"$name [UUID(${workflowId.shortString})$jobTag]"

      override def slf4jLoggers: Set[Logger] = Set.empty
    }

    override lazy val callEngineFunctions = functions
  }

  private val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)
  private val workingDisk = JesWorkingDisk(DiskType.SSD, 200)

  val DockerAndDiskRuntime =
    """
      |runtime {
      |  docker: "ubuntu:latest"
      |  disks: "local-disk 200 SSD"
      |}
    """.stripMargin

  private def buildPreemptibleJobDescriptor(attempt: Int, preemptible: Int): BackendJobDescriptor = {
    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      NamespaceWithWorkflow.load(YoSup.replace("[PREEMPTIBLE]", s"preemptible: $preemptible")),
      Inputs,
      NoOptions
    )

    val job = workflowDescriptor.workflowNamespace.workflow.calls.head
    val key = BackendJobDescriptorKey(job, None, attempt)
    val runtimeAttributes = makeRuntimeAttributes(job)
    BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Inputs)
  }

  private def executionActor(jobDescriptor: BackendJobDescriptor,
                             configurationDescriptor: BackendConfigurationDescriptor,
                             promise: Promise[BackendJobExecutionResponse],
                             errorCode: Int,
                             innerErrorCode: Int): ActorRef = {

    // Mock/stub out the bits that would reach out to JES.
    val run = mock[Run]
    run.status() returns Failed(errorCode, Option(s"$innerErrorCode: I seen some things man"), Seq.empty, Option("fakeMachine"), Option("fakeZone"), Option("fakeInstance"))

    val handle = JesPendingExecutionHandle(jobDescriptor, Seq.empty, run, None)

    class ExecuteOrRecoverActor extends TestableJesJobExecutionActor(jobDescriptor, promise, jesConfiguration) {
      override def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[ExecutionHandle] = Future.successful(handle)
    }

    system.actorOf(Props(new ExecuteOrRecoverActor), "ExecuteOrRecoverActor-" + UUID.randomUUID)
  }

  private def run(attempt: Int, preemptible: Int, errorCode: Int, innerErrorCode: Int): BackendJobExecutionResponse = {
    within(Timeout) {
      val promise = Promise[BackendJobExecutionResponse]()
      val jobDescriptor =  buildPreemptibleJobDescriptor(attempt, preemptible)
      val backend = executionActor(jobDescriptor, JesBackendConfigurationDescriptor, promise, errorCode, innerErrorCode)
      backend ! Execute
      Await.result(promise.future, Timeout)
    }
  }

  def buildPreemptibleTestActorRef(attempt: Int, preemptible: Int): TestActorRef[TestableJesJobExecutionActor] = {
    val jobDescriptor = buildPreemptibleJobDescriptor(attempt, preemptible)
    val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
    TestActorRef(props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")
  }

  behavior of "JesAsyncBackendJobExecutionActor"

  { // Set of "handle call failures appropriately with respect to preemption" tests
    val expectations = Table(
      ("attempt", "preemptible", "errorCode", "innerErrorCode", "shouldRetry"),
      // No preemptible attempts allowed, nothing should be retryable.
      (1, 0, 10, 13, false),
      (1, 0, 10, 14, false),
      (1, 0, 10, 15, false),
      (1, 0, 11, 13, false),
      (1, 0, 11, 14, false),
      // 1 preemptible attempt allowed, but not all failures represent preemptions.
      (1, 1, 10, 13, true),
      (1, 1, 10, 14, true),
      (1, 1, 10, 15, false),
      (1, 1, 11, 13, false),
      (1, 1, 11, 14, false),
      // 1 preemptible attempt allowed, but now on the second attempt nothing should be retryable.
      (2, 1, 10, 13, false),
      (2, 1, 10, 14, false),
      (2, 1, 10, 15, false),
      (2, 1, 11, 13, false),
      (2, 1, 11, 14, false)
    )

    expectations foreach { case (attempt, preemptible, errorCode, innerErrorCode, shouldRetry) =>
      it should s"handle call failures appropriately with respect to preemption (attempt=$attempt, preemptible=$preemptible, errorCode=$errorCode, innerErrorCode=$innerErrorCode)" in {
        run(attempt, preemptible, errorCode, innerErrorCode).getClass.getSimpleName match {
          case "FailedNonRetryableResponse" => false shouldBe shouldRetry
          case "FailedRetryableResponse" => true shouldBe shouldRetry
          case huh => fail(s"Unexpected response class name: '$huh'")
        }
      }
    }
  }

  it should "not restart 2 of 1 unexpected shutdowns without another preemptible VM" in {
    val actorRef = buildPreemptibleTestActorRef(2, 1)
    val jesBackend = actorRef.underlyingActor
    val handle = mock[JesPendingExecutionHandle]
    implicit val ec = system.dispatcher

    val failedStatus = Failed(10, Some("14: VM XXX shut down unexpectedly."), Seq.empty, Option("fakeMachine"), Option("fakeZone"), Option("fakeInstance"))
    val executionResult = Await.result(jesBackend.executionResult(failedStatus, handle), 2.seconds)
    executionResult.isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    val failedHandle = executionResult.asInstanceOf[FailedNonRetryableExecutionHandle]
    failedHandle.returnCode shouldBe None
  }

  it should "restart 1 of 1 unexpected shutdowns without another preemptible VM" in {
    val actorRef = buildPreemptibleTestActorRef(1, 1)
    val jesBackend = actorRef.underlyingActor
    val handle = mock[JesPendingExecutionHandle]
    implicit val ec = system.dispatcher

    val failedStatus = Failed(10, Some("14: VM XXX shut down unexpectedly."), Seq.empty, Option("fakeMachine"), Option("fakeZone"), Option("fakeInstance"))
    val executionResult = Await.result(jesBackend.executionResult(failedStatus, handle), 2.seconds)
    executionResult.isInstanceOf[FailedRetryableExecutionHandle] shouldBe true
    val retryableHandle = executionResult.asInstanceOf[FailedRetryableExecutionHandle]
    retryableHandle.throwable.isInstanceOf[PreemptedException] shouldBe true
    retryableHandle.returnCode shouldBe None
    val preemptedException = retryableHandle.throwable.asInstanceOf[PreemptedException]
    preemptedException.getMessage should include("will be restarted with a non-preemptible VM")
  }

  it should "restart 1 of 2 unexpected shutdowns with another preemptible VM" in {
    val actorRef = buildPreemptibleTestActorRef(1, 2)
    val jesBackend = actorRef.underlyingActor
    val handle = mock[JesPendingExecutionHandle]
    implicit val ec = system.dispatcher

    val failedStatus = Failed(10, Some("14: VM XXX shut down unexpectedly."), Seq.empty, Option("fakeMachine"), Option("fakeZone"), Option("fakeInstance"))
    val executionResult = Await.result(jesBackend.executionResult(failedStatus, handle), 2.seconds)
    executionResult.isInstanceOf[FailedRetryableExecutionHandle] shouldBe true
    val retryableHandle = executionResult.asInstanceOf[FailedRetryableExecutionHandle]
    retryableHandle.throwable.isInstanceOf[PreemptedException] shouldBe true
    retryableHandle.returnCode shouldBe None
    val preemptedException2 = retryableHandle.throwable.asInstanceOf[PreemptedException]
    preemptedException2.getMessage should include("will be restarted with another preemptible VM")
  }

  it should "handle Failure Status for various errors" in {
    val actorRef = buildPreemptibleTestActorRef(1, 1)
    val jesBackend = actorRef.underlyingActor
    val handle = mock[JesPendingExecutionHandle]
    implicit val ec = system.dispatcher

    Await.result(jesBackend.executionResult(
      Failed(10, Some("15: Other type of error."), Seq.empty, Option("fakeMachine"), Option("fakeZone"), Option("fakeInstance")), handle), 2.seconds
    ).isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    Await.result(jesBackend.executionResult(
      Failed(11, Some("14: Wrong errorCode."), Seq.empty, Option("fakeMachine"), Option("fakeZone"), Option("fakeInstance")), handle), 2.seconds
    ).isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    Await.result(jesBackend.executionResult(
      Failed(10, Some("Weird error message."), Seq.empty, Option("fakeMachine"), Option("fakeZone"), Option("fakeInstance")), handle), 2.seconds
    ).isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    Await.result(jesBackend.executionResult(
      Failed(10, Some("UnparsableInt: Even weirder error message."), Seq.empty, Option("fakeMachine"), Option("fakeZone"), Option("fakeInstance")), handle), 2.seconds
    ).isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    Await.result(jesBackend.executionResult(
      Failed(10, None, Seq.empty, Option("fakeMachine"), Option("fakeZone"), Option("fakeInstance")), handle), 2.seconds
    ).isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    Await.result(jesBackend.executionResult(
      Failed(10, Some("Operation canceled at"), Seq.empty, Option("fakeMachine"), Option("fakeZone"), Option("fakeInstance")), handle), 2.seconds
    ) shouldBe AbortedExecutionHandle

    actorRef.stop()
  }

  it should "map GCS paths and *only* GCS paths to local" in {
    val stringKey = "abc"
    val stringVal = WdlString("abc")
    val localFileKey = "lf"
    val localFileVal = WdlFile("/blah/abc")
    val gcsFileKey = "gcsf"
    val gcsFileVal = WdlFile("gs://blah/abc")

    val inputs = Map(
      stringKey -> stringVal,
      localFileKey -> localFileVal,
      gcsFileKey -> gcsFileVal
    )

    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      NamespaceWithWorkflow.load(YoSup.replace("[PREEMPTIBLE]", "")),
      inputs,
      NoOptions
    )

    val key = BackendJobDescriptorKey(workflowDescriptor.workflowNamespace.workflow.calls.head, None, 1)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, Map.empty, inputs)

    val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
    val testActorRef = TestActorRef[TestableJesJobExecutionActor](
      props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

    val mappedInputs = jobDescriptor.inputs mapValues testActorRef.underlyingActor.gcsPathToLocal

    mappedInputs(stringKey) match {
      case WdlString(v) => assert(v.equalsIgnoreCase(stringVal.value))
      case _ => fail("test setup error")
    }

    mappedInputs(localFileKey) match {
      case wdlFile: WdlFile => assert(wdlFile.value.equalsIgnoreCase(localFileVal.value))
      case _ => fail("test setup error")
    }

    mappedInputs(gcsFileKey) match {
      case wdlFile: WdlFile => assert(wdlFile.value.equalsIgnoreCase("blah/abc"))
      case _ => fail("test setup error")
    }
  }

  it should "generate correct JesFileInputs from a WdlMap" in {
    val inputs = Map(
      "stringToFileMap" -> WdlMap(WdlMapType(WdlStringType, WdlFileType), Map(
        WdlString("stringTofile1") -> WdlFile("gs://path/to/stringTofile1"),
        WdlString("stringTofile2") -> WdlFile("gs://path/to/stringTofile2")
      )),
      "fileToStringMap" -> WdlMap(WdlMapType(WdlFileType, WdlStringType), Map(
        WdlFile("gs://path/to/fileToString1") -> WdlString("fileToString1"),
        WdlFile("gs://path/to/fileToString2") -> WdlString("fileToString2")
      )),
      "fileToFileMap" -> WdlMap(WdlMapType(WdlFileType, WdlFileType), Map(
        WdlFile("gs://path/to/fileToFile1Key") -> WdlFile("gs://path/to/fileToFile1Value"),
        WdlFile("gs://path/to/fileToFile2Key") -> WdlFile("gs://path/to/fileToFile2Value")
      )),
      "stringToString" -> WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(
        WdlString("stringToString1") -> WdlString("path/to/stringToString1"),
        WdlString("stringToString2") -> WdlString("path/to/stringToString2")
      ))
    )

    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      NamespaceWithWorkflow.load(SampleWdl.CurrentDirectory.asWorkflowSources(DockerAndDiskRuntime).wdlSource),
      inputs,
      NoOptions
    )

    val job = workflowDescriptor.workflowNamespace.workflow.calls.head
    val runtimeAttributes = makeRuntimeAttributes(job)
    val key = BackendJobDescriptorKey(job, None, 1)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, inputs)

    val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
    val testActorRef = TestActorRef[TestableJesJobExecutionActor](
      props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

    val jesInputs = testActorRef.underlyingActor.generateJesInputs(jobDescriptor)
    jesInputs should have size 8
    jesInputs should contain(JesFileInput(
      "stringToFileMap-0", "gs://path/to/stringTofile1", Paths.get("path/to/stringTofile1"), workingDisk))
    jesInputs should contain(JesFileInput(
      "stringToFileMap-1", "gs://path/to/stringTofile2", Paths.get("path/to/stringTofile2"), workingDisk))
    jesInputs should contain(JesFileInput(
      "fileToStringMap-0", "gs://path/to/fileToString1", Paths.get("path/to/fileToString1"), workingDisk))
    jesInputs should contain(JesFileInput(
      "fileToStringMap-1", "gs://path/to/fileToString2", Paths.get("path/to/fileToString2"), workingDisk))
    jesInputs should contain(JesFileInput(
      "fileToFileMap-0", "gs://path/to/fileToFile1Key", Paths.get("path/to/fileToFile1Key"), workingDisk))
    jesInputs should contain(JesFileInput(
      "fileToFileMap-1", "gs://path/to/fileToFile1Value", Paths.get("path/to/fileToFile1Value"), workingDisk))
    jesInputs should contain(JesFileInput(
      "fileToFileMap-2", "gs://path/to/fileToFile2Key", Paths.get("path/to/fileToFile2Key"), workingDisk))
    jesInputs should contain(JesFileInput(
      "fileToFileMap-3", "gs://path/to/fileToFile2Value", Paths.get("path/to/fileToFile2Value"), workingDisk))
  }

  def makeJesActorRef(sampleWdl: SampleWdl, callName: LocallyQualifiedName, inputs: Map[FullyQualifiedName, WdlValue],
                      functions: JesExpressionFunctions = TestableJesExpressionFunctions):
  TestActorRef[TestableJesJobExecutionActor] = {
    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      NamespaceWithWorkflow.load(sampleWdl.asWorkflowSources(DockerAndDiskRuntime).wdlSource),
      inputs,
      NoOptions
    )

    val call = workflowDescriptor.workflowNamespace.workflow.findCallByName(callName).get
    val key = BackendJobDescriptorKey(call, None, 1)
    val runtimeAttributes = makeRuntimeAttributes(call)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, inputs)

    val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration, functions))
    TestActorRef[TestableJesJobExecutionActor](props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")
  }

  it should "generate correct JesOutputs" in {
    val inputs = Map(
      "in" -> WdlFile("gs://a/b/c.txt")
    )
    val jesBackend = makeJesActorRef(SampleWdl.FilePassingWorkflow, "a", inputs).underlyingActor
    val jobDescriptor = jesBackend.jobDescriptor
    val workflowId = jesBackend.workflowId
    val jesInputs = jesBackend.generateJesInputs(jobDescriptor)
    jesInputs should have size 1
    jesInputs should contain(JesFileInput("in-0", "gs://a/b/c.txt", Paths.get("a/b/c.txt"), workingDisk))
    val jesOutputs = jesBackend.generateJesOutputs(jobDescriptor)
    jesOutputs should have size 1
    jesOutputs should contain(JesFileOutput("out",
      s"gs://my-cromwell-workflows-bucket/file_passing/$workflowId/call-a/out", Paths.get("out"), workingDisk))
  }

  it should "generate correct JesInputs when a command line contains a write_lines call in it" in {
    val inputs = Map(
      "strs" -> WdlArray(WdlArrayType(WdlStringType), Seq("A", "B", "C").map(WdlString))
    )

    class TestJesExpressionFunctions extends JesExpressionFunctions(
      List(MockGcsFileSystemBuilder.mockGcsFileSystem), TestableCallContext) {
      override def write_lines(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
        Success(WdlFile(s"gs://some/path/file.txt"))
      }
    }

    val functions = new TestJesExpressionFunctions
    val jesBackend = makeJesActorRef(SampleWdl.ArrayIO, "serialize", inputs, functions).underlyingActor
    val jobDescriptor = jesBackend.jobDescriptor
    val jesInputs = jesBackend.generateJesInputs(jobDescriptor)
    jesInputs should have size 1
    jesInputs should contain(JesFileInput(
      "c6fd5c91-0", "gs://some/path/file.txt", Paths.get("some/path/file.txt"), workingDisk))
    val jesOutputs = jesBackend.generateJesOutputs(jobDescriptor)
    jesOutputs should have size 0
  }

  it should "generate correct JesFileInputs from a WdlArray" in {
    val inputs = Map(
      "fileArray" ->
        WdlArray(WdlArrayType(WdlFileType), Seq(WdlFile("gs://path/to/file1"), WdlFile("gs://path/to/file2")))
    )

    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      NamespaceWithWorkflow.load(SampleWdl.CurrentDirectory.asWorkflowSources(DockerAndDiskRuntime).wdlSource),
      inputs,
      NoOptions
    )

    val job = workflowDescriptor.workflowNamespace.workflow.calls.head
    val runtimeAttributes = makeRuntimeAttributes(job)
    val key = BackendJobDescriptorKey(job, None, 1)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, inputs)

    val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
    val testActorRef = TestActorRef[TestableJesJobExecutionActor](
      props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

    val jesInputs = testActorRef.underlyingActor.generateJesInputs(jobDescriptor)
    jesInputs should have size 2
    jesInputs should contain(JesFileInput("fileArray-0", "gs://path/to/file1", Paths.get("path/to/file1"), workingDisk))
    jesInputs should contain(JesFileInput("fileArray-1", "gs://path/to/file2", Paths.get("path/to/file2"), workingDisk))
  }

  it should "generate correct JesFileInputs from a WdlFile" in {
    val inputs = Map(
      "file1" -> WdlFile("gs://path/to/file1"),
      "file2" -> WdlFile("gs://path/to/file2")
    )

    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      NamespaceWithWorkflow.load(SampleWdl.CurrentDirectory.asWorkflowSources(DockerAndDiskRuntime).wdlSource),
      inputs,
      NoOptions
    )

    val job = workflowDescriptor.workflowNamespace.workflow.calls.head
    val runtimeAttributes = makeRuntimeAttributes(job)
    val key = BackendJobDescriptorKey(job, None, 1)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, inputs)

    val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
    val testActorRef = TestActorRef[TestableJesJobExecutionActor](
      props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

    val jesInputs = testActorRef.underlyingActor.generateJesInputs(jobDescriptor)
    jesInputs should have size 2
    jesInputs should contain(JesFileInput("file1-0", "gs://path/to/file1", Paths.get("path/to/file1"), workingDisk))
    jesInputs should contain(JesFileInput("file2-0", "gs://path/to/file2", Paths.get("path/to/file2"), workingDisk))
  }

  it should "convert local Paths back to corresponding GCS paths in JesOutputs" in {
    val jesOutputs = Seq(
      JesFileOutput("/cromwell_root/path/to/file1", "gs://path/to/file1",
        Paths.get("/cromwell_root/path/to/file1"), workingDisk),
      JesFileOutput("/cromwell_root/path/to/file2", "gs://path/to/file2",
        Paths.get("/cromwell_root/path/to/file2"), workingDisk),
      JesFileOutput("/cromwell_root/path/to/file3", "gs://path/to/file3",
        Paths.get("/cromwell_root/path/to/file3"), workingDisk),
      JesFileOutput("/cromwell_root/path/to/file4", "gs://path/to/file4",
        Paths.get("/cromwell_root/path/to/file4"), workingDisk),
      JesFileOutput("/cromwell_root/path/to/file5", "gs://path/to/file5",
        Paths.get("/cromwell_root/path/to/file5"), workingDisk)
    )
    val outputValues = Seq(
      WdlFile("/cromwell_root/path/to/file1"),
      WdlArray(WdlArrayType(WdlFileType), Seq(
        WdlFile("/cromwell_root/path/to/file2"), WdlFile("/cromwell_root/path/to/file3"))),
      WdlMap(WdlMapType(WdlFileType, WdlFileType), Map(
        WdlFile("/cromwell_root/path/to/file4") -> WdlFile("/cromwell_root/path/to/file5")
      ))
    )

    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      NamespaceWithWorkflow.load(SampleWdl.EmptyString.asWorkflowSources(DockerAndDiskRuntime).wdlSource),
      Map.empty,
      NoOptions
    )

    val key = BackendJobDescriptorKey(workflowDescriptor.workflowNamespace.workflow.calls.head, None, 1)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, Map.empty, Map.empty)

    val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
    val testActorRef = TestActorRef[TestableJesJobExecutionActor](
      props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

    val result = outputValues map testActorRef.underlyingActor.wdlValueToGcsPath(jesOutputs)
    result should have size 3
    result should contain(WdlFile("gs://path/to/file1"))
    result should contain(WdlArray(WdlArrayType(WdlFileType),
      Seq(WdlFile("gs://path/to/file2"), WdlFile("gs://path/to/file3")))
    )
    result should contain(WdlMap(WdlMapType(WdlFileType, WdlFileType),
      Map(WdlFile("gs://path/to/file4") -> WdlFile("gs://path/to/file5")))
    )
  }

  it should "create a JesFileInput for the monitoring script, when specified" in {
    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      NamespaceWithWorkflow.load(SampleWdl.EmptyString.asWorkflowSources(DockerAndDiskRuntime).wdlSource),
      Map.empty,
      WorkflowOptions.fromJsonString("""{"monitoring_script": "gs://path/to/script"}""").get
    )

    val job = workflowDescriptor.workflowNamespace.workflow.calls.head
    val runtimeAttributes = makeRuntimeAttributes(job)
    val key = BackendJobDescriptorKey(job, None, 1)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty)

    val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
    val testActorRef = TestActorRef[TestableJesJobExecutionActor](
      props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

    testActorRef.underlyingActor.monitoringScript shouldBe
      Some(JesFileInput("monitoring-in", "gs://path/to/script", Paths.get("monitoring.sh"), workingDisk))
  }

  it should "not create a JesFileInput for the monitoring script, when not specified" in {
    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      NamespaceWithWorkflow.load(SampleWdl.EmptyString.asWorkflowSources(DockerAndDiskRuntime).wdlSource),
      Map.empty,
      NoOptions
    )

    val key = BackendJobDescriptorKey(workflowDescriptor.workflowNamespace.workflow.calls.head, None, 1)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, Map.empty, Map.empty)

    val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
    val testActorRef = TestActorRef[TestableJesJobExecutionActor](
      props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

    testActorRef.underlyingActor.monitoringScript shouldBe None
  }

  it should "return JES log paths for non-scattered call" in {
    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId(UUID.fromString("e6236763-c518-41d0-9688-432549a8bf7c")),
      NamespaceWithWorkflow.load(
        SampleWdl.HelloWorld.asWorkflowSources(""" runtime {docker: "ubuntu:latest"} """).wdlSource),
      Map.empty,
      WorkflowOptions.fromJsonString(""" {"jes_gcs_root": "gs://path/to/gcs_root"} """).get
    )

    val call = workflowDescriptor.workflowNamespace.workflow.findCallByName("hello").get
    val key = BackendJobDescriptorKey(call, None, 1)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, Map.empty, Map.empty)

    val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
    val testActorRef = TestActorRef[TestableJesJobExecutionActor](
      props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

    val jesBackend = testActorRef.underlyingActor

    // TODO: NioGcsPath.equals not implemented, so use toString instead
    jesBackend.jesCallPaths.stdoutPath should be(a[NioGcsPath])
    jesBackend.jesCallPaths.stdoutPath.toString shouldBe
      "gs://path/to/gcs_root/hello/e6236763-c518-41d0-9688-432549a8bf7c/call-hello/hello-stdout.log"
    jesBackend.jesCallPaths.stderrPath should be(a[NioGcsPath])
    jesBackend.jesCallPaths.stderrPath.toString shouldBe
      "gs://path/to/gcs_root/hello/e6236763-c518-41d0-9688-432549a8bf7c/call-hello/hello-stderr.log"
    jesBackend.jesCallPaths.jesLogPath should be(a[NioGcsPath])
    jesBackend.jesCallPaths.jesLogPath.toString shouldBe
      "gs://path/to/gcs_root/hello/e6236763-c518-41d0-9688-432549a8bf7c/call-hello/hello.log"
  }

  it should "return JES log paths for scattered call" in {
    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId(UUID.fromString("e6236763-c518-41d0-9688-432549a8bf7d")),
      NamespaceWithWorkflow.load(
        new SampleWdl.ScatterWdl().asWorkflowSources(""" runtime {docker: "ubuntu:latest"} """).wdlSource),
      Map.empty,
      WorkflowOptions.fromJsonString(""" {"jes_gcs_root": "gs://path/to/gcs_root"} """).get
    )

    val call = workflowDescriptor.workflowNamespace.workflow.findCallByName("B").get
    val key = BackendJobDescriptorKey(call, Option(2), 1)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, Map.empty, Map.empty)

    val props = Props(new TestableJesJobExecutionActor(jobDescriptor, Promise(), jesConfiguration))
    val testActorRef = TestActorRef[TestableJesJobExecutionActor](
      props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

    val jesBackend = testActorRef.underlyingActor

    jesBackend.jesCallPaths.stdoutPath should be(a[NioGcsPath])
    jesBackend.jesCallPaths.stdoutPath.toString shouldBe
      "gs://path/to/gcs_root/w/e6236763-c518-41d0-9688-432549a8bf7d/call-B/shard-2/B-2-stdout.log"
    jesBackend.jesCallPaths.stderrPath should be(a[NioGcsPath])
    jesBackend.jesCallPaths.stderrPath.toString shouldBe
      "gs://path/to/gcs_root/w/e6236763-c518-41d0-9688-432549a8bf7d/call-B/shard-2/B-2-stderr.log"
    jesBackend.jesCallPaths.jesLogPath should be(a[NioGcsPath])
    jesBackend.jesCallPaths.jesLogPath.toString shouldBe
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

  private def makeRuntimeAttributes(job: Call) = {
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(job.task.runtimeAttributes, TestableJesExpressionFunctions, Map.empty)
    RuntimeAttributeDefinition.addDefaultsToAttributes(JesBackendLifecycleActorFactory.staticRuntimeAttributeDefinitions, NoOptions)(evaluatedAttributes.get) // Fine to throw the exception if this "get" fails. This is a test after all!
  }
}
