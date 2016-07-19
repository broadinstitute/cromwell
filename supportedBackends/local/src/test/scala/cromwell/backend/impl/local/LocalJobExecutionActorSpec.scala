package cromwell.backend.impl.local

import java.nio.file.{Files, Paths}

import akka.testkit.TestActorRef
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, FailedNonRetryableResponse, SucceededResponse}
import cromwell.backend.io.TestWorkflows._
import cromwell.backend.io.{WorkflowPaths, JobPaths, TestWorkflows}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobDescriptorKey, BackendSpec}
import cromwell.core._
import org.scalatest.FlatSpecLike
import org.scalatest.concurrent.PatienceConfiguration._
import org.scalatest.mock.MockitoSugar
import org.scalatest.prop.TableDrivenPropertyChecks
import wdl4s.types._
import wdl4s.values._

import scala.concurrent.duration._

class LocalJobExecutionActorSpec extends TestKitSuite("LocalJobExecutionActorSpec") with FlatSpecLike
  with BackendSpec with MockitoSugar with TestFileUtil with TableDrivenPropertyChecks {

  val globalConfig = ConfigFactory.load()
  val backendConfig = globalConfig.getConfig("backend.providers.Local.config")
  val defaultBackendConfigDescriptor = new BackendConfigurationDescriptor(backendConfig, globalConfig)

  def localBackend(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor) = {
    val workflowPaths = new WorkflowPaths(jobDescriptor.descriptor, configurationDescriptor.backendConfig)
    TestActorRef(new LocalJobExecutionActor(jobDescriptor, configurationDescriptor, LocalBackendInitializationData(workflowPaths), scala.concurrent.ExecutionContext.global)).underlyingActor
  }

  behavior of "LocalBackend"

  it should "execute an hello world workflow" in {
    val expectedOutputs: JobOutputs = Map(
      "salutation" -> JobOutput(WdlString("Hello you !"), None)
    )
    val expectedResponse = SucceededResponse(mock[BackendJobDescriptorKey], Some(0), expectedOutputs)
    val wf = new TestWorkflow(buildWorkflowDescriptor(HelloWorld), defaultBackendConfigDescriptor, expectedResponse)
    val backend = localBackend(jobDescriptorFromSingleCallWorkflow(wf.workflowDescriptor), wf.config)
    testWorkflow(wf, backend)
  }

  it should "execute an hello world workflow on Docker" taggedAs DockerTest in {
    val expectedOutputs: JobOutputs = Map(
      "salutation" -> JobOutput(WdlString("Hello you !"), None)
    )
    val expectedResponse = SucceededResponse(mock[BackendJobDescriptorKey], Some(0), expectedOutputs)
    val wf = new TestWorkflow(buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" }"""), defaultBackendConfigDescriptor, expectedResponse)
    val backend = localBackend(jobDescriptorFromSingleCallWorkflow(wf.workflowDescriptor), wf.config)
    testWorkflow(wf, backend)
  }

  it should "send back an execution failure if the task fails" in {
    val expectedResponse = FailedNonRetryableResponse(mock[BackendJobDescriptorKey], new Exception(""), Option(1))
    val wf = new TestWorkflow(buildWorkflowDescriptor(GoodbyeWorld), defaultBackendConfigDescriptor, expectedResponse)
    val backend = localBackend(jobDescriptorFromSingleCallWorkflow(wf.workflowDescriptor), wf.config)
    testWorkflow(wf, backend)
  }

  def localizationSpec(docker: Boolean) = {
    def templateConf(localizers: String) = BackendConfigurationDescriptor(
      ConfigFactory.parseString(
        s"""{
            |  root = "local-cromwell-executions"
            |  filesystems {
            |    local {
            |      localization = [
            |        $localizers
            |      ]
            |    }
            |  }
            |}
        """.stripMargin),
      ConfigFactory.parseString("{}")
    )

    val hardConf = templateConf("hard-link")
    val symConf = templateConf("soft-link")
    val copyConf = templateConf("copy")

    val jsonInputFile = createCannedFile("localize", "content from json inputs").toPath.toAbsolutePath.toString
    val callInputFile = createCannedFile("localize", "content from call inputs").toPath.toAbsolutePath.toString
    val inputs = Map(
      "inputFileFromCallInputs" -> WdlFile(callInputFile),
      "inputFileFromJson" -> WdlFile(jsonInputFile)
    )

    val expectedOutputs: JobOutputs = Map(
      "out" -> JobOutput(WdlArray(WdlArrayType(WdlStringType),
        Array(
          WdlString("content from json inputs"),
          WdlString("content from call inputs"))
      ), None)
    )

    val confs = List(
      (hardConf, false),
      (copyConf, false)
    ) ++ (if (!docker) List((symConf, true)) else List.empty)

    val localizers = Table(
      ("conf", "isSymLink"),
      confs:_*
    )

    forAll(localizers) { (conf, isSymlink) =>
      val runtime = if (docker) """runtime { docker: "ubuntu:latest" } """ else ""
      val wf = buildWorkflowDescriptor(InputFiles, inputs, runtime = runtime)
      val backend = localBackend(jobDescriptorFromSingleCallWorkflow(wf, inputs), conf)
      val jobDescriptor: BackendJobDescriptor = jobDescriptorFromSingleCallWorkflow(wf)
      val expectedResponse = SucceededResponse(jobDescriptor.key, Some(0), expectedOutputs)

      val jobPaths = new JobPaths(wf, conf.backendConfig, jobDescriptor.key)

      whenReady(backend.execute) { executionResponse =>
        assertResponse(executionResponse, expectedResponse)
        val localizedJsonInputFile = Paths.get(jobPaths.callRoot.toString, jsonInputFile.toString)
        val localizedCallInputFile = Paths.get(jobPaths.callRoot.toString, callInputFile.toString)

        Files.isSymbolicLink(localizedJsonInputFile) shouldBe isSymlink
        val realJsonInputFile = if (isSymlink) Files.readSymbolicLink(localizedJsonInputFile).toFile else localizedJsonInputFile.toFile
        realJsonInputFile.exists() shouldBe true

        Files.isSymbolicLink(localizedCallInputFile) shouldBe isSymlink
        val realCallInputFile = if (isSymlink) Files.readSymbolicLink(localizedJsonInputFile).toFile else localizedCallInputFile.toFile
        realCallInputFile.exists() shouldBe true
      }
    }
  }

  it should "execute calls with input files and localize them appropriately" in {
    localizationSpec(docker = false)
  }

  it should "execute calls with input files and localize them appropriately (in Docker)" taggedAs DockerTest in {
    localizationSpec(docker = true)
  }

  it should "abort a job and kill a process" in {
    val wf = buildWorkflowDescriptor(Sleep10)
    val jobDescriptor: BackendJobDescriptor = jobDescriptorFromSingleCallWorkflow(wf)
    val backend = localBackend(jobDescriptor, defaultBackendConfigDescriptor)

    val execute = backend.execute
    // TODO: PBE: This test needs work. If the abort fires to quickly, it causes a race condition in waitAndPostProcess.
    Thread.sleep(1000L)
    //backend.process shouldNot be(empty) <-- Currently cannot access the private var backend.process
    backend.abort

    // TODO: PBE: abort doesn't actually seem to abort. It runs the full 10 seconsds, then returns the response.
    whenReady(execute, Timeout(10.seconds)) { executionResponse =>
      executionResponse shouldBe a[AbortedResponse]
    }
  }

  it should "execute shards from a scatter" in {
    val wf = buildWorkflowDescriptor(TestWorkflows.Scatter)

    val call = wf.workflowNamespace.workflow.calls.head

    0 to 2 foreach { shard =>
      // This assumes that engine will give us the evaluated value of the scatter item at the correct index
      // If this is not the case, more context/logic will need to be moved to the backend so it can figure it out by itself
      val symbolMaps: Map[LocallyQualifiedName, WdlInteger] = Map("intNumber" -> WdlInteger(shard))

      val jd: BackendJobDescriptor = new BackendJobDescriptor(wf, new BackendJobDescriptorKey(call, Option(shard), 1), symbolMaps)
      val backend = localBackend(jd, defaultBackendConfigDescriptor)
      val response = SucceededResponse(mock[BackendJobDescriptorKey], Some(0), Map("out" -> JobOutput(WdlInteger(shard), None)))
      executeJobAndAssertOutputs(backend, response)
    }
  }

  it should "post process outputs" in {
    val inputFile = createCannedFile("localize", "content from json inputs").toPath.toAbsolutePath.toString
    val inputs = Map {
      "inputFile" -> WdlFile(inputFile)
    }
    val wf = buildWorkflowDescriptor(OutputProcess, inputs)
    val jobDescriptor: BackendJobDescriptor = jobDescriptorFromSingleCallWorkflow(wf, inputs)
    val backend = localBackend(jobDescriptor, defaultBackendConfigDescriptor)
    val jobPaths = new JobPaths(wf, defaultBackendConfigDescriptor.backendConfig, jobDescriptor.key)
    val expectedA = WdlFile(jobPaths.callRoot.resolve("a").toAbsolutePath.toString)
    val expectedB = WdlFile(jobPaths.callRoot.resolve("dir").toAbsolutePath.resolve("b").toString)
    val expectedOutputs = Map (
      "o1" -> JobOutput(expectedA, None),
      "o2" -> JobOutput(
        WdlArray(WdlArrayType(WdlFileType), Seq(expectedA, expectedB)), None
      ),
      "o3" -> JobOutput(WdlFile(inputFile), None)
    )
    val expectedResponse = SucceededResponse(jobDescriptor.key, Some(0), expectedOutputs)

    executeJobAndAssertOutputs(backend, expectedResponse)
  }

  it should "fail post processing if an output fail is not found" in {
    val expectedResponse = FailedNonRetryableResponse(mock[BackendJobDescriptorKey], new Throwable("Failed post processing of outputs"), Option(0))
    val wf = new TestWorkflow(buildWorkflowDescriptor(MissingOutputProcess), defaultBackendConfigDescriptor, expectedResponse)
    val backend = localBackend(jobDescriptorFromSingleCallWorkflow(wf.workflowDescriptor), wf.config)
    testWorkflow(wf, backend)
  }
}
