package cromwell.backend.sfs

import java.nio.file.{Files, Path, Paths}

import akka.actor.Props
import akka.testkit.TestActorRef
import better.files._
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, BackendJobExecutionResponse, FailedNonRetryableResponse, SucceededResponse}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.io.TestWorkflows._
import cromwell.backend.io.{JobPaths, TestWorkflows, WorkflowPaths}
import cromwell.backend.validation.{DockerValidation, RuntimeAttributesValidation}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobDescriptorKey, BackendSpec}
import cromwell.core.Tags._
import cromwell.core._
import org.scalatest.FlatSpecLike
import org.scalatest.concurrent.PatienceConfiguration._
import org.scalatest.mock.MockitoSugar
import org.scalatest.prop.TableDrivenPropertyChecks
import wdl4s.types._
import wdl4s.util.AggregatedException
import wdl4s.values._

import scala.concurrent.Promise
import scala.concurrent.duration._

class TestLocalAsyncJobExecutionActor(override val params: SharedFileSystemAsyncJobExecutionActorParams)
  extends BackgroundAsyncJobExecutionActor {
  override lazy val processArgs = {
    val script = jobPaths.script.fullPath
    if (isDockerRun) {
      val docker = RuntimeAttributesValidation.extract(DockerValidation.instance, validatedRuntimeAttributes)
      val cwd = jobPaths.callRoot.fullPath
      val dockerCwd = jobPaths.callDockerRoot.fullPath
      SharedFileSystemCommand("/bin/bash", "-c",
        s"docker run --rm -v $cwd:$dockerCwd -i $docker /bin/bash < $script")
    } else {
      SharedFileSystemCommand("/bin/bash", script)
    }
  }
}

class SharedFileSystemJobExecutionActorSpec extends TestKitSuite("SharedFileSystemJobExecutionActorSpec")
  with FlatSpecLike with BackendSpec with MockitoSugar with TableDrivenPropertyChecks {

  def createBackend(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor) = {
    createBackendRef(jobDescriptor, configurationDescriptor).underlyingActor
  }

  def createBackendRef(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor) = {
    val workflowPaths = new WorkflowPaths(jobDescriptor.workflowDescriptor, configurationDescriptor.backendConfig)
    val initializationData = new SharedFileSystemBackendInitializationData(workflowPaths,
      SharedFileSystemValidatedRuntimeAttributesBuilder.default.withValidation(DockerValidation.optional))

    def propsCreator(completionPromise: Promise[BackendJobExecutionResponse]): Props = {
      val params = SharedFileSystemAsyncJobExecutionActorParams(emptyActor, jobDescriptor,
        configurationDescriptor, completionPromise, Option(initializationData))
      Props(classOf[TestLocalAsyncJobExecutionActor], params)
    }

    TestActorRef(new SharedFileSystemJobExecutionActor(
      jobDescriptor, configurationDescriptor, propsCreator))
  }

  behavior of "SharedFileSystemJobExecutionActor"

  it should "execute an hello world workflow" in {
    val expectedOutputs: JobOutputs = Map(
      "salutation" -> JobOutput(WdlString("Hello you !"))
    )
    val expectedResponse = SucceededResponse(mock[BackendJobDescriptorKey], Some(0), expectedOutputs)
    val workflow = TestWorkflow(buildWorkflowDescriptor(HelloWorld), emptyBackendConfig, expectedResponse)
    val backend = createBackend(jobDescriptorFromSingleCallWorkflow(workflow.workflowDescriptor), workflow.config)
    testWorkflow(workflow, backend)
  }

  it should "execute an hello world workflow on Docker" taggedAs DockerTest in {
    val expectedOutputs: JobOutputs = Map(
      "salutation" -> JobOutput(WdlString("Hello you !"))
    )
    val expectedResponse = SucceededResponse(mock[BackendJobDescriptorKey], Some(0), expectedOutputs)
    val workflow = TestWorkflow(buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" }"""),
      emptyBackendConfig, expectedResponse)
    val backend = createBackend(jobDescriptorFromSingleCallWorkflow(workflow.workflowDescriptor), workflow.config)
    testWorkflow(workflow, backend)
  }

  it should "send back an execution failure if the task fails" in {
    val expectedResponse = FailedNonRetryableResponse(mock[BackendJobDescriptorKey], new Exception(""), Option(1))
    val workflow = TestWorkflow(buildWorkflowDescriptor(GoodbyeWorld), emptyBackendConfig, expectedResponse)
    val backend = createBackend(jobDescriptorFromSingleCallWorkflow(workflow.workflowDescriptor), workflow.config)
    testWorkflow(workflow, backend)
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

    val jsonInputFile = createCannedFile("localize", "content from json inputs").fullPath
    val callInputFile = createCannedFile("localize", "content from call inputs").fullPath
    val inputs = Map(
      "inputFileFromCallInputs" -> WdlFile(callInputFile),
      "inputFileFromJson" -> WdlFile(jsonInputFile)
    )

    val expectedOutputs: JobOutputs = Map(
      "out" -> JobOutput(WdlArray(WdlArrayType(WdlStringType),
        Array(
          WdlString("content from json inputs"),
          WdlString("content from call inputs")))))

    val confs = List(
      (hardConf, false),
      (copyConf, false)
    ) ++ (if (!docker) List((symConf, true)) else List.empty)

    val localizers = Table(
      ("conf", "isSymLink"),
      confs: _*
    )

    forAll(localizers) { (conf, isSymlink) =>
      val runtime = if (docker) """runtime { docker: "ubuntu:latest" } """ else ""
      val workflowDescriptor = buildWorkflowDescriptor(InputFiles, inputs, runtime = runtime)
      val backend = createBackend(jobDescriptorFromSingleCallWorkflow(workflowDescriptor, inputs), conf)
      val jobDescriptor: BackendJobDescriptor = jobDescriptorFromSingleCallWorkflow(workflowDescriptor)
      val expectedResponse = SucceededResponse(jobDescriptor.key, Some(0), expectedOutputs)

      val jobPaths = new JobPaths(workflowDescriptor, conf.backendConfig, jobDescriptor.key)

      whenReady(backend.execute) { executionResponse =>
        assertResponse(executionResponse, expectedResponse)
        val localizedJsonInputFile = Paths.get(jobPaths.callRoot.toString, jsonInputFile.toString)
        val localizedCallInputFile = Paths.get(jobPaths.callRoot.toString, callInputFile.toString)

        Files.isSymbolicLink(localizedJsonInputFile) shouldBe isSymlink
        val realJsonInputFile =
          if (isSymlink) Files.readSymbolicLink(localizedJsonInputFile) else localizedJsonInputFile
        realJsonInputFile.toFile should exist

        Files.isSymbolicLink(localizedCallInputFile) shouldBe isSymlink
        val realCallInputFile =
          if (isSymlink) Files.readSymbolicLink(localizedJsonInputFile) else localizedCallInputFile
        realCallInputFile.toFile should exist
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
    val workflowDescriptor = buildWorkflowDescriptor(Sleep10)
    val jobDescriptor: BackendJobDescriptor = jobDescriptorFromSingleCallWorkflow(workflowDescriptor)
    val backendRef = createBackendRef(jobDescriptor, emptyBackendConfig)
    val backend = backendRef.underlyingActor

    val execute = backend.execute
    // TODO: PBE: This test needs work. If the abort fires to quickly, it causes a race condition in waitAndPostProcess.
    Thread.sleep(1000L)
    //backend.process shouldNot be(empty) <-- Currently cannot access the private var backend.process
    backendRef ! AbortJobCommand

    // TODO: PBE: abort doesn't actually seem to abort. It runs the full 10 seconsds, then returns the response.
    whenReady(execute, Timeout(10.seconds)) { executionResponse =>
      executionResponse shouldBe a[AbortedResponse]
    }
  }

  it should "execute shards from a scatter" in {
    val workflowDescriptor = buildWorkflowDescriptor(TestWorkflows.Scatter)

    val call = workflowDescriptor.workflowNamespace.workflow.calls.head

    0 to 2 foreach { shard =>
      // This assumes that engine will give us the evaluated value of the scatter item at the correct index
      // If this is not the case, more context/logic will need to be moved to the backend so it can figure it out by itself
      val symbolMaps: Map[LocallyQualifiedName, WdlInteger] = Map("intNumber" -> WdlInteger(shard))

      val jobDescriptor: BackendJobDescriptor =
        BackendJobDescriptor(workflowDescriptor, BackendJobDescriptorKey(call, Option(shard), 1), Map.empty, symbolMaps)
      val backend = createBackend(jobDescriptor, emptyBackendConfig)
      val response =
        SucceededResponse(mock[BackendJobDescriptorKey], Some(0), Map("out" -> JobOutput(WdlInteger(shard))))
      executeJobAndAssertOutputs(backend, response)
    }
  }

  it should "post process outputs" in {
    val inputFile = createCannedFile("localize", "content from json inputs").fullPath
    val inputs = Map {
      "inputFile" -> WdlFile(inputFile)
    }
    val workflowDescriptor = buildWorkflowDescriptor(OutputProcess, inputs)
    val jobDescriptor: BackendJobDescriptor = jobDescriptorFromSingleCallWorkflow(workflowDescriptor, inputs)
    val backend = createBackend(jobDescriptor, emptyBackendConfig)
    val jobPaths = new JobPaths(workflowDescriptor, emptyBackendConfig.backendConfig, jobDescriptor.key)
    val expectedA = WdlFile(jobPaths.callRoot.resolve("a").toAbsolutePath.toString)
    val expectedB = WdlFile(jobPaths.callRoot.resolve("dir").toAbsolutePath.resolve("b").toString)
    val expectedOutputs = Map(
      "o1" -> JobOutput(expectedA),
      "o2" -> JobOutput(WdlArray(WdlArrayType(WdlFileType), Seq(expectedA, expectedB))),
      "o3" -> JobOutput(WdlFile(inputFile))
    )
    val expectedResponse = SucceededResponse(jobDescriptor.key, Some(0), expectedOutputs)

    executeJobAndAssertOutputs(backend, expectedResponse)
  }

  it should "fail post processing if an output fail is not found" in {
    val expectedResponse = FailedNonRetryableResponse(mock[BackendJobDescriptorKey],
      AggregatedException(Seq.empty, "Could not process output, file not found"), Option(0))
    val workflow = TestWorkflow(buildWorkflowDescriptor(MissingOutputProcess), emptyBackendConfig, expectedResponse)
    val backend = createBackend(jobDescriptorFromSingleCallWorkflow(workflow.workflowDescriptor), workflow.config)
    testWorkflow(workflow, backend)
  }

  def createCannedFile(prefix: String, contents: String, dir: Option[Path] = None): Path = {
    val suffix = ".out"
    val file = dir match {
      case Some(path) => Files.createTempFile(path, prefix, suffix)
      case None => Files.createTempFile(prefix, suffix)
    }
    file.write(contents)
    file
  }
}
