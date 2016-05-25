package cromwell.backend.impl.local

import java.nio.file.{Files, Paths}

import akka.testkit.TestActorRef
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, FailedNonRetryableResponse, SucceededResponse}
import cromwell.backend.io.{TestWorkflows, JobPaths, BackendTestkitSpec}
import cromwell.backend.io.BackendTestkitSpec._
import cromwell.backend.io.TestWorkflows._
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobDescriptorKey}
import cromwell.core._
import org.scalatest.FlatSpec
import org.scalatest.mock.MockitoSugar
import org.scalatest.prop.TableDrivenPropertyChecks
import wdl4s.types._
import wdl4s.values._

class LocalJobExecutionActorSpec extends FlatSpec with BackendTestkitSpec with MockitoSugar with TestFileUtil with TableDrivenPropertyChecks {

  val globalConfig = ConfigFactory.load()
  val backendConfig = globalConfig.getConfig("backend.providers.Local.config")
  val defaultBackendConfigDescriptor = new BackendConfigurationDescriptor(backendConfig, globalConfig)

  def localBackend(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor) =
    TestActorRef(new LocalJobExecutionActor(jobDescriptor, configurationDescriptor)).underlyingActor

  behavior of "LocalBackend"

  it should "execute an hello world workflow" in {
    val expectedOutputs: CallOutputs = Map(
      "salutation" -> CallOutput(WdlString("Hello you !"), None)
    )
    val expectedResponse = SucceededResponse(mock[BackendJobDescriptorKey], expectedOutputs)
    val wf = new TestWorkflow(buildWorkflowDescriptor(HelloWorld), defaultBackendConfigDescriptor, expectedResponse)
    val backend = localBackend(jobDescriptorFromSingleCallWorkflow(wf.workflowDescriptor), wf.config)
    testWorkflow(wf, backend)
  }

  it should "execute an hello world workflow on Docker" taggedAs DockerTest in {
    val expectedOutputs: CallOutputs = Map(
      "salutation" -> CallOutput(WdlString("Hello you !"), None)
    )
    val expectedResponse = SucceededResponse(mock[BackendJobDescriptorKey], expectedOutputs)
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

  it should "execute calls with input files and localize them appropriately" in {

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

    val expectedOutputs: CallOutputs = Map(
      "out" -> CallOutput(WdlArray(WdlArrayType(WdlStringType),
        Array(
          WdlString("content from json inputs"),
          WdlString("content from call inputs"))
      ), None)
    )

    val localizers = Table(
      ("conf", "isSymLink"),
      (hardConf, false),
      (copyConf, false),
      (symConf, true)
    )

    forAll(localizers) { (conf, isSymlink) =>
      val wf = buildWorkflowDescriptor(InputFiles, inputs)
      val backend = localBackend(jobDescriptorFromSingleCallWorkflow(wf, inputs), conf)
      val jobDescriptor: BackendJobDescriptor = jobDescriptorFromSingleCallWorkflow(wf)
      val expectedResponse = SucceededResponse(jobDescriptor.key, expectedOutputs)

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

  it should "abort a job and kill a process" in {
    val wf = buildWorkflowDescriptor(Sleep10)
    val jobDescriptor: BackendJobDescriptor = jobDescriptorFromSingleCallWorkflow(wf)
    val backend = localBackend(jobDescriptor, defaultBackendConfigDescriptor)

    val execute = backend.execute
    val abort = backend.abortJob

    whenReady(execute) { executionResponse =>
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
      val response = SucceededResponse(mock[BackendJobDescriptorKey], Map("out" -> CallOutput(WdlInteger(shard), None)))
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
      "o1" -> CallOutput(expectedA, None),
      "o2" -> CallOutput(
        WdlArray(WdlArrayType(WdlFileType), Seq(expectedA, expectedB)), None
      ),
      "o3" -> CallOutput(WdlFile(inputFile), None)
    )
    val expectedResponse = SucceededResponse(jobDescriptor.key, expectedOutputs)

    executeJobAndAssertOutputs(backend, expectedResponse)
  }

  it should "fail post processing if an output fail is not found" in {
    val expectedResponse = FailedNonRetryableResponse(mock[BackendJobDescriptorKey], new Throwable("Failed post processing of outputs"), Option(0))
    val wf = new TestWorkflow(buildWorkflowDescriptor(MissingOutputProcess), defaultBackendConfigDescriptor, expectedResponse)
    val backend = localBackend(jobDescriptorFromSingleCallWorkflow(wf.workflowDescriptor), wf.config)
    testWorkflow(wf, backend)
  }
}
