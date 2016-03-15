package cromwell.engine.backend.jes

import java.net.URL
import java.nio.file.{FileSystems, Paths}
import java.util.UUID

import com.google.api.client.testing.http.{HttpTesting, MockHttpTransport, MockLowLevelHttpRequest, MockLowLevelHttpResponse}
import cromwell.CromwellTestkitSpec
import cromwell.engine._
import cromwell.engine.backend.io.filesystem.gcs.{GcsFileSystem, NioGcsPath}
import cromwell.engine.backend.jes.JesBackend.{JesFileInput, JesFileOutput}
import cromwell.engine.backend.jes.Run.Failed
import cromwell.engine.backend.jes.authentication._
import cromwell.engine.backend.runtimeattributes.{CromwellRuntimeAttributes, DiskType}
import cromwell.engine.backend.{AbortedExecutionHandle, BackendCallJobDescriptor, FailedExecutionHandle, RetryableExecutionHandle}
import cromwell.engine.io.gcs._
import cromwell.engine.workflow.{BackendCallKey, WorkflowOptions}
import cromwell.util.{EncryptionSpec, SampleWdl}
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import org.slf4j.{Logger, LoggerFactory}
import org.specs2.mock.Mockito
import wdl4s.types.{WdlArrayType, WdlFileType, WdlMapType, WdlStringType}
import wdl4s.values._
import wdl4s.{Call, CallInputs, Task}

import scala.concurrent.Await
import scala.concurrent.duration._
import scala.util.{Success, Try}

class JesBackendSpec extends FlatSpec with Matchers with Mockito with BeforeAndAfterAll {
  val testWorkflowManagerSystem = new CromwellTestkitSpec.TestWorkflowManagerSystem()
  val actorSystem = testWorkflowManagerSystem.actorSystem
  val workingDisk = JesWorkingDisk(DiskType.SSD, 200)

  override protected def afterAll() = {
    testWorkflowManagerSystem.shutdownTestActorSystem()
    super.afterAll()
  }

  val clientSecrets = SimpleClientSecrets("id", "secrets")
  val jesBackend = new JesBackend(actorSystem) {
    private val anyString = ""
    private val anyURL: URL = null
    override lazy val jesConf = new JesAttributes(
      project = anyString,
      executionBucket = anyString,
      endpointUrl = anyURL,
      maxPollingInterval = 600) {
    }
    override def jesUserConnection(workflow: WorkflowDescriptor) = null
    override lazy val jesCromwellInterface = null
    override lazy val googleConf = GoogleConfiguration("appName", ServiceAccountMode("accountID", "pem"), Option(Refresh(clientSecrets)))
  }

  "executionResult" should "handle Failure Status" in {
    import scala.concurrent.ExecutionContext.Implicits.global

    val logger: Logger = LoggerFactory.getLogger("JesBackendSpecLogger")
    val wd = mock[WorkflowDescriptor]
    wd.workflowLogger returns logger
    val task = mock[Task]
    task.outputs returns Seq.empty
    val call = mock[Call]
    call.task returns task

    class PreemptionJobDescriptor(attempt: Int, max: Int) extends BackendCallJobDescriptor(wd, BackendCallKey(call, None, attempt)) {
      private val attributes = mock[CromwellRuntimeAttributes].preemptible returns max
      override def callRuntimeAttributes: CromwellRuntimeAttributes = attributes
    }

    val handle = mock[JesPendingExecutionHandle]
    handle.jobDescriptor returns new PreemptionJobDescriptor(attempt = 2, max = 1)

    val executionResult0 = Await.result(jesBackend.executionResult(new Failed(10, Some("14: VM XXX shut down unexpectedly."), Seq.empty), handle), 2.seconds)
    executionResult0.isInstanceOf[FailedExecutionHandle] shouldBe true
    val failedHandle0 = executionResult0.asInstanceOf[FailedExecutionHandle]
    failedHandle0.returnCode shouldBe None

    handle.jobDescriptor returns new PreemptionJobDescriptor(attempt = 1, max = 1)

    val executionResult = Await.result(jesBackend.executionResult(new Failed(10, Some("14: VM XXX shut down unexpectedly."), Seq.empty), handle), 2.seconds)
    executionResult.isInstanceOf[RetryableExecutionHandle] shouldBe true
    val retryableHandle = executionResult.asInstanceOf[RetryableExecutionHandle]
    retryableHandle.throwable.isInstanceOf[PreemptedException] shouldBe true
    retryableHandle.returnCode shouldBe None
    val preemptedException = retryableHandle.throwable.asInstanceOf[PreemptedException]
    preemptedException.getMessage should include ("will be restarted with a non-preemptible VM")

    handle.jobDescriptor returns new PreemptionJobDescriptor(attempt = 1, max = 2)

    val executionResult2 = Await.result(jesBackend.executionResult(new Failed(10, Some("14: VM XXX shut down unexpectedly."), Seq.empty), handle), 2.seconds)
    executionResult2.isInstanceOf[RetryableExecutionHandle] shouldBe true
    val retryableHandle2 = executionResult2.asInstanceOf[RetryableExecutionHandle]
    retryableHandle2.throwable.isInstanceOf[PreemptedException] shouldBe true
    retryableHandle2.returnCode shouldBe None
    val preemptedException2 = retryableHandle2.throwable.asInstanceOf[PreemptedException]
    preemptedException2.getMessage should include ("will be restarted with another preemptible VM")

    Await.result(jesBackend.executionResult(new Failed(10, Some("15: Other type of error."), Seq.empty), handle), 2.seconds).isInstanceOf[FailedExecutionHandle] shouldBe true
    Await.result(jesBackend.executionResult(new Failed(11, Some("14: Wrong errorCode."), Seq.empty), handle), 2.seconds).isInstanceOf[FailedExecutionHandle] shouldBe true
    Await.result(jesBackend.executionResult(new Failed(10, Some("Weird error message."), Seq.empty), handle), 2.seconds).isInstanceOf[FailedExecutionHandle] shouldBe true
    Await.result(jesBackend.executionResult(new Failed(10, Some("UnparsableInt: Even weirder error message."), Seq.empty), handle), 2.seconds).isInstanceOf[FailedExecutionHandle] shouldBe true
    Await.result(jesBackend.executionResult(new Failed(10, None, Seq.empty), handle), 2.seconds).isInstanceOf[FailedExecutionHandle] shouldBe true
    Await.result(jesBackend.executionResult(new Failed(10, Some("Operation canceled at"), Seq.empty), handle), 2.seconds) shouldBe AbortedExecutionHandle
  }

  "JesBackend" should "consider 403 as a fatal exception" in {
    val transport = new MockHttpTransport() {
      override def buildRequest(method: String, url: String) = {
        new MockLowLevelHttpRequest() {
          override def execute() = new MockLowLevelHttpResponse().setStatusCode(403)
        }
      }
    }
    val request = transport.createRequestFactory().buildGetRequest(HttpTesting.SIMPLE_GENERIC_URL)
    val mockedResponse = Try(request.execute()).failed.get
    JesBackend.isFatalJesException(mockedResponse) shouldBe true
  }

  it should "consider 429 as a transient exception" in {
    val transport = new MockHttpTransport() {
      override def buildRequest(method: String, url: String) = {
        new MockLowLevelHttpRequest() {
          override def execute() = new MockLowLevelHttpResponse().setStatusCode(429)
        }
      }
    }
    val request = transport.createRequestFactory().buildGetRequest(HttpTesting.SIMPLE_GENERIC_URL)
    val mockedResponse = Try(request.execute()).failed.get
    JesBackend.isTransientJesException(mockedResponse) shouldBe true
  }

  "adjustInputPaths" should "map GCS paths and *only* GCS paths to local" in {
    val ignoredCall = mock[BackendCallKey]
    val stringKey = "abc"
    val stringVal = WdlString("abc")
    val localFileKey = "lf"
    val localFileVal = WdlFile("/blah/abc")
    val gcsFileKey = "gcsf"
    val gcsFileVal = WdlFile("gs://blah/abc")
    val emptyRuntimeAttributes = CromwellRuntimeAttributes.defaults

    val inputs: CallInputs = collection.immutable.HashMap(
      stringKey -> stringVal,
      localFileKey -> localFileVal,
      gcsFileKey -> gcsFileVal
    )

    // This should only ever be used in this test to grab some locallyQualifiedInputs. So leave the rest null:
    val jobDescriptor = BackendCallJobDescriptor(null, null, inputs)

    val mockedBackendCall = mock[JesBackendCall]
    mockedBackendCall.jobDescriptor returns jobDescriptor
    val mappedInputs: CallInputs = new JesBackend(actorSystem).adjustInputPaths(mockedBackendCall.jobDescriptor)

    mappedInputs.get(stringKey).get match {
      case WdlString(v) => assert(v.equalsIgnoreCase(stringVal.value))
      case _ => fail("test setup error")
    }

    mappedInputs.get(localFileKey).get match {
      case wdlFile: WdlFile => assert(wdlFile.value.equalsIgnoreCase(localFileVal.value))
      case _ => fail("test setup error")
    }

    mappedInputs.get(gcsFileKey).get match {
      case wdlFile: WdlFile => assert(wdlFile.value.equalsIgnoreCase("blah/abc"))
      case _ => fail("test setup error")
    }
  }

  "workflow options existence" should "be verified when localizing with Refresh Token" in {
    EncryptionSpec.assumeAes256Cbc()

    val goodOptions = WorkflowOptions.fromMap(Map("refresh_token" -> "token")).get
    val missingToken = WorkflowOptions.fromMap(Map.empty).get

    try {
      jesBackend.assertWorkflowOptions(goodOptions)
    } catch {
      case e: IllegalArgumentException => fail("Correct options validation should not throw an exception.")
      case t: Throwable =>
        t.printStackTrace()
        fail(s"Unexpected exception: ${t.getMessage}")
    }

    the [IllegalArgumentException] thrownBy {
      jesBackend.assertWorkflowOptions(missingToken)
    } should have message s"Missing parameters in workflow options: refresh_token"
  }

  it should "create a GcsAuthInformation instance" in {
    val workflowDescriptor = mock[WorkflowDescriptor]
    val mockedWfOptions = mock[WorkflowOptions]
    workflowDescriptor.workflowOptions returns mockedWfOptions
    mockedWfOptions.get("refresh_token") returns Success("myRefreshToken")

    jesBackend.getGcsAuthInformation(workflowDescriptor) shouldBe Some(GcsLocalizing(clientSecrets, "myRefreshToken"))
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
    val descriptor = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), SampleWdl.CurrentDirectory.asWorkflowSources())
    val jobDescriptor = mock[BackendCallJobDescriptor]
    jobDescriptor.locallyQualifiedInputs returns inputs
    val runtimeAttributes = mock[CromwellRuntimeAttributes]
    runtimeAttributes.disks returns Seq(workingDisk)
    jobDescriptor.callRuntimeAttributes returns runtimeAttributes
    val key = BackendCallKey(descriptor.namespace.workflow.findCallByName("whereami").get, None, 1)
    jobDescriptor.key returns key
    val jesInputs = jesBackend.generateJesInputs(jobDescriptor)
    jesInputs should have size 8
    jesInputs should contain(JesFileInput("stringToFileMap-0", "gs://path/to/stringTofile1", Paths.get("path/to/stringTofile1"), workingDisk))
    jesInputs should contain(JesFileInput("stringToFileMap-1", "gs://path/to/stringTofile2", Paths.get("path/to/stringTofile2"), workingDisk))
    jesInputs should contain(JesFileInput("fileToStringMap-0", "gs://path/to/fileToString1", Paths.get("path/to/fileToString1"), workingDisk))
    jesInputs should contain(JesFileInput("fileToStringMap-1", "gs://path/to/fileToString2", Paths.get("path/to/fileToString2"), workingDisk))
    jesInputs should contain(JesFileInput("fileToFileMap-0", "gs://path/to/fileToFile1Key", Paths.get("path/to/fileToFile1Key"), workingDisk))
    jesInputs should contain(JesFileInput("fileToFileMap-1", "gs://path/to/fileToFile1Value", Paths.get("path/to/fileToFile1Value"), workingDisk))
    jesInputs should contain(JesFileInput("fileToFileMap-2", "gs://path/to/fileToFile2Key", Paths.get("path/to/fileToFile2Key"), workingDisk))
    jesInputs should contain(JesFileInput("fileToFileMap-3", "gs://path/to/fileToFile2Value", Paths.get("path/to/fileToFile2Value"), workingDisk))
  }

  def makeJobDescriptor(wdl: SampleWdl,
                        callName: String,
                        inputs: Map[String, WdlValue],
                        lookup: String => WdlValue,
                        functions: JesCallEngineFunctions = new JesCallEngineFunctions(List(GcsFileSystem.defaultGcsFileSystem), new CallContext("root", "out", "err"))): BackendCallJobDescriptor = {

    val descriptor = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), wdl.asWorkflowSources()).copy(wfContext = new WorkflowContext("gs://foobar"))
    val jobDescriptor = mock[BackendCallJobDescriptor]
    val runtimeAttributes = mock[CromwellRuntimeAttributes]
    runtimeAttributes.disks returns Seq(workingDisk)
    jobDescriptor.workflowDescriptor returns descriptor
    jobDescriptor.key returns new BackendCallKey(descriptor.namespace.workflow.findCallByName(callName).get, None, 1)
    jobDescriptor.locallyQualifiedInputs returns inputs

    val call = descriptor.namespace.workflow.findCallByName(callName).get
    jobDescriptor.key returns BackendCallKey(call, None, 1)
    // This will be less ugly once BackendCall is really deleted and the Filesystems PR is merged
    jobDescriptor.callRootPath returns new NioGcsPath(Seq("call", "gcs", "path").toArray, absolute=true, isDirectory = true)(GcsFileSystem.defaultGcsFileSystem)
    jobDescriptor.lookupFunction(inputs) returns lookup
    jobDescriptor.callRuntimeAttributes returns runtimeAttributes
    jobDescriptor.callEngineFunctions returns functions
    jobDescriptor
  }

  it should "generate correct JesOutputs" in {
    val inputs = Map(
      "in" -> WdlFile("gs://a/b/c.txt")
    )
    val jobDescriptor = makeJobDescriptor(SampleWdl.FilePassingWorkflow, "a", inputs, (s: String) => inputs.get(s).get)
    val jesInputs = jesBackend.generateJesInputs(jobDescriptor)
    jesInputs should have size 1
    jesInputs should contain(JesFileInput("in-0", "gs://a/b/c.txt", Paths.get("a/b/c.txt"), workingDisk))
    val jesOutputs = jesBackend.generateJesOutputs(jobDescriptor)
    jesOutputs should have size 1
    jesOutputs should contain(JesFileOutput("out", "gs://call/gcs/path/out", Paths.get("out"), workingDisk))
  }

  it should "generate correct JesInputs when a command line contains a write_lines call in it" in {
    val inputs = Map(
      "strs" -> WdlArray(WdlArrayType(WdlStringType), Seq("A", "B", "C").map(WdlString))
    )
    class TestEngineFunctions(context: CallContext) extends JesCallEngineFunctions(List(GcsFileSystem.defaultGcsFileSystem), context) {
      override def write_lines(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
        Success(WdlFile(s"gs://some/path/file.txt"))
      }
    }
    val functions = new TestEngineFunctions(new CallContext("root", "stdout", "stderr"))
    val jobDescriptor = makeJobDescriptor(SampleWdl.ArrayIO, "serialize", inputs, (s: String) => inputs.get(s).get, functions)
    val jesInputs = jesBackend.generateJesInputs(jobDescriptor)
    jesInputs should have size 1
    jesInputs should contain(JesFileInput("c6fd5c91-0", "gs://some/path/file.txt", Paths.get("some/path/file.txt"), workingDisk))
    val jesOutputs = jesBackend.generateJesOutputs(jobDescriptor)
    jesOutputs should have size 0
  }

  it should "generate correct JesFileInputs from a WdlArray" in {
    val inputs = Map(
      "fileArray" -> WdlArray(WdlArrayType(WdlFileType), Seq(WdlFile("gs://path/to/file1"), WdlFile("gs://path/to/file2")))
    )
    val descriptor = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), SampleWdl.CurrentDirectory.asWorkflowSources())
    val jobDescriptor = mock[BackendCallJobDescriptor]
    jobDescriptor.locallyQualifiedInputs returns inputs
    val runtimeAttributes = mock[CromwellRuntimeAttributes]
    runtimeAttributes.disks returns Seq(workingDisk)
    jobDescriptor.callRuntimeAttributes returns runtimeAttributes
    val key = BackendCallKey(descriptor.namespace.workflow.findCallByName("whereami").get, None, 1)
    jobDescriptor.key returns key

    val jesInputs = jesBackend.generateJesInputs(jobDescriptor)
    jesInputs should have size 2
    jesInputs should contain(JesFileInput("fileArray-0", "gs://path/to/file1", Paths.get("path/to/file1"), workingDisk))
    jesInputs should contain(JesFileInput("fileArray-1", "gs://path/to/file2", Paths.get("path/to/file2"), workingDisk))
  }

  it should "generate correct JesFileInputs from a WdlFile" in {
    val inputs = Map(
      "file1" -> WdlFile("gs://path/to/file1"),
      "file2" -> WdlFile("gs://path/to/file2")
    )
    val descriptor = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), SampleWdl.CurrentDirectory.asWorkflowSources())
    val jobDescriptor = mock[BackendCallJobDescriptor]
    jobDescriptor.locallyQualifiedInputs returns inputs
    val runtimeAttributes = mock[CromwellRuntimeAttributes]
    runtimeAttributes.disks returns Seq(workingDisk)
    jobDescriptor.callRuntimeAttributes returns runtimeAttributes
    val key = BackendCallKey(descriptor.namespace.workflow.findCallByName("whereami").get, None, 1)
    jobDescriptor.key returns key

    val jesInputs = jesBackend.generateJesInputs(jobDescriptor)
    jesInputs should have size 2
    jesInputs should contain(JesFileInput("file1-0", "gs://path/to/file1", Paths.get("path/to/file1"), workingDisk))
    jesInputs should contain(JesFileInput("file2-0", "gs://path/to/file2", Paths.get("path/to/file2"), workingDisk))
  }

  it should "convert local Paths back to corresponding GCS paths in JesOutputs" in {
    val jesOutputs = Seq(
      JesFileOutput("/cromwell_root/path/to/file1", "gs://path/to/file1", Paths.get("/cromwell_root/path/to/file1"), workingDisk),
      JesFileOutput("/cromwell_root/path/to/file2", "gs://path/to/file2", Paths.get("/cromwell_root/path/to/file2"), workingDisk),
      JesFileOutput("/cromwell_root/path/to/file3", "gs://path/to/file3", Paths.get("/cromwell_root/path/to/file3"), workingDisk),
      JesFileOutput("/cromwell_root/path/to/file4", "gs://path/to/file4", Paths.get("/cromwell_root/path/to/file4"), workingDisk),
      JesFileOutput("/cromwell_root/path/to/file5", "gs://path/to/file5", Paths.get("/cromwell_root/path/to/file5"), workingDisk)
    )
    val outputValues = Seq(
      WdlFile("/cromwell_root/path/to/file1"),
      WdlArray(WdlArrayType(WdlFileType), Seq(WdlFile("/cromwell_root/path/to/file2"), WdlFile("/cromwell_root/path/to/file3"))),
      WdlMap(WdlMapType(WdlFileType, WdlFileType), Map(
        WdlFile("/cromwell_root/path/to/file4") -> WdlFile("/cromwell_root/path/to/file5")
      ))
    )
    val result = outputValues map jesBackend.wdlValueToGcsPath(jesOutputs)
    result should have size 3
    result should contain(WdlFile("gs://path/to/file1"))
    result should contain(WdlArray(WdlArrayType(WdlFileType),
      Seq(WdlFile("gs://path/to/file2"), WdlFile("gs://path/to/file3")))
    )
    result should contain(WdlMap(WdlMapType(WdlFileType, WdlFileType),
      Map(WdlFile("gs://path/to/file4") -> WdlFile("gs://path/to/file5")))
    )
  }

  it should "create a JesFileInput for the monitoring script, if specified" in {
    val jobDescriptor = mock[BackendCallJobDescriptor]
    val wd = mock[WorkflowDescriptor]
    jobDescriptor.workflowDescriptor returns wd
    val runtimeAttributes = mock[CromwellRuntimeAttributes]
    runtimeAttributes.disks returns Seq(workingDisk)
    jobDescriptor.callRuntimeAttributes returns runtimeAttributes

    wd.workflowOptions returns WorkflowOptions.fromJsonString("""{"monitoring_script": "gs://path/to/script"}""").get
    JesBackend.monitoringIO(jobDescriptor) shouldBe Some(JesFileInput("monitoring-in", "gs://path/to/script", Paths.get("monitoring.sh"), workingDisk))

    wd.workflowOptions returns WorkflowOptions.fromJsonString("""{}""").get
    JesBackend.monitoringIO(jobDescriptor) shouldBe None
  }

  "JesBackendCall" should "return JES log paths for non-scattered call" in {
    val wd = WorkflowDescriptor(WorkflowId(UUID.fromString("e6236763-c518-41d0-9688-432549a8bf7c")), SampleWdl.HelloWorld.asWorkflowSources(
      runtime = """ runtime {docker: "ubuntu:latest"} """,
      workflowOptions = """ {"jes_gcs_root": "gs://path/to/gcs_root"} """
    )).copy(wfContext = new WorkflowContext("gs://path/to/gcs_root")).copy(
      fileSystems = List(GcsFileSystem.defaultGcsFileSystem, FileSystems.getDefault)
    )

    val call = wd.namespace.workflow.findCallByName("hello").get
    val backendCall = jesBackend.bindCall(BackendCallJobDescriptor(wd, BackendCallKey(call, None, 1)))
    val stdoutstderr = backendCall.stdoutStderr

    stdoutstderr.stdout shouldBe WdlFile("gs://path/to/gcs_root/hello/e6236763-c518-41d0-9688-432549a8bf7c/call-hello/hello-stdout.log")
    stdoutstderr.stderr shouldBe WdlFile("gs://path/to/gcs_root/hello/e6236763-c518-41d0-9688-432549a8bf7c/call-hello/hello-stderr.log")

    stdoutstderr.backendLogs shouldBe defined
    val logsMap = stdoutstderr.backendLogs.get
    logsMap should contain key "log"
    logsMap("log") shouldBe WdlFile("gs://path/to/gcs_root/hello/e6236763-c518-41d0-9688-432549a8bf7c/call-hello/hello.log")
  }

  it should "return JES log paths for scattered call" in {
    val wd = WorkflowDescriptor(WorkflowId(UUID.fromString("e6236763-c518-41d0-9688-432549a8bf7c")), new SampleWdl.ScatterWdl().asWorkflowSources(
      runtime = """ runtime {docker: "ubuntu:latest"} """,
      workflowOptions = """ {"jes_gcs_root": "gs://path/to/gcs_root"} """
    )).copy(wfContext = new WorkflowContext("gs://path/to/gcs_root")).copy(
      fileSystems = List(GcsFileSystem.defaultGcsFileSystem, FileSystems.getDefault)
    )
    val call = wd.namespace.workflow.findCallByName("B").get
    val backendCall = jesBackend.bindCall(BackendCallJobDescriptor(wd, BackendCallKey(call, Some(2), 1)))
    val stdoutstderr = backendCall.stdoutStderr

    stdoutstderr.stdout shouldBe WdlFile("gs://path/to/gcs_root/w/e6236763-c518-41d0-9688-432549a8bf7c/call-B/shard-2/B-2-stdout.log")
    stdoutstderr.stderr shouldBe WdlFile("gs://path/to/gcs_root/w/e6236763-c518-41d0-9688-432549a8bf7c/call-B/shard-2/B-2-stderr.log")

    stdoutstderr.backendLogs shouldBe defined
    val logsMap = stdoutstderr.backendLogs.get
    logsMap should contain key "log"
    logsMap("log") shouldBe WdlFile("gs://path/to/gcs_root/w/e6236763-c518-41d0-9688-432549a8bf7c/call-B/shard-2/B-2.log")
  }
}
