package cromwell.backend.sfs

import java.nio.file.{Files, Paths}

import akka.testkit.TestDuration
import better.files._
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, JobFailedNonRetryableResponse, JobSucceededResponse}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.io.TestWorkflows._
import cromwell.backend.io.{JobPathsWithDocker, TestWorkflows}
import cromwell.backend.sfs.TestLocalAsyncJobExecutionActor._
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobDescriptorKey, BackendSpec, RuntimeAttributeDefinition}
import cromwell.core.Tags._
import cromwell.core._
import cromwell.services.keyvalue.KeyValueServiceActor.{KvJobKey, KvPair, ScopedKey}
import lenthall.exception.AggregatedException
import org.scalatest.concurrent.PatienceConfiguration.Timeout
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{Assertion, FlatSpecLike, OptionValues}
import wdl4s.types._
import wdl4s.values._

import scala.concurrent.duration._

class SharedFileSystemJobExecutionActorSpec extends TestKitSuite("SharedFileSystemJobExecutionActorSpec")
  with FlatSpecLike with BackendSpec with TableDrivenPropertyChecks with OptionValues {

  behavior of "SharedFileSystemJobExecutionActor"

  lazy val runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition] =
    SharedFileSystemValidatedRuntimeAttributesBuilder.default.definitions.toSet

  def executeSpec(docker: Boolean): Any = {
    val expectedOutputs: CallOutputs = Map(
      "salutation" -> JobOutput(WdlString("Hello you !"))
    )
    val expectedResponse = JobSucceededResponse(mock[BackendJobDescriptorKey], Some(0), expectedOutputs, None, Seq.empty)
    val runtime = if (docker) """runtime { docker: "ubuntu:latest" }""" else ""
    val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = runtime)
    val workflow = TestWorkflow(workflowDescriptor, emptyBackendConfig, expectedResponse)
    val backend = createBackend(jobDescriptorFromSingleCallWorkflow(workflow.workflowDescriptor, Map.empty, WorkflowOptions.empty, runtimeAttributeDefinitions), workflow.config)
    testWorkflow(workflow, backend)
  }

  it should "execute an hello world workflow" in {
    executeSpec(docker = false)
  }

  it should "execute an hello world workflow on Docker" taggedAs DockerTest in {
    executeSpec(docker = true)
  }

  it should "send back an execution failure if the task fails" in {
    val expectedResponse =
      JobFailedNonRetryableResponse(mock[BackendJobDescriptorKey], new RuntimeException(""), Option(1))
    val workflow = TestWorkflow(buildWorkflowDescriptor(GoodbyeWorld), emptyBackendConfig, expectedResponse)
    val backend = createBackend(jobDescriptorFromSingleCallWorkflow(workflow.workflowDescriptor, Map.empty, WorkflowOptions.empty, runtimeAttributeDefinitions), workflow.config)
    testWorkflow(workflow, backend)
  }

  def localizationSpec(docker: Boolean): Assertion = {
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

    val jsonInputFile = createCannedFile("localize", "content from json inputs").pathAsString
    val callInputFile = createCannedFile("localize", "content from call inputs").pathAsString
    val inputs = Map(
      "wf_localize.localize.inputFileFromCallInputs" -> WdlFile(callInputFile),
      "wf_localize.localize.inputFileFromJson" -> WdlFile(jsonInputFile)
    )

    val expectedOutputs: CallOutputs = Map(
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
      val backend = createBackend(jobDescriptorFromSingleCallWorkflow(workflowDescriptor, inputs, WorkflowOptions.empty, runtimeAttributeDefinitions), conf)
      val jobDescriptor: BackendJobDescriptor = jobDescriptorFromSingleCallWorkflow(workflowDescriptor, inputs, WorkflowOptions.empty, runtimeAttributeDefinitions)
      val expectedResponse = JobSucceededResponse(jobDescriptor.key, Some(0), expectedOutputs, None, Seq.empty)

      val jobPaths = new JobPathsWithDocker(jobDescriptor.key, workflowDescriptor, conf.backendConfig)

      whenReady(backend.execute) { executionResponse =>
        assertResponse(executionResponse, expectedResponse)
        val localizedJsonInputFile = Paths.get(jobPaths.callInputsRoot.toString, jsonInputFile)
        val localizedCallInputFile = Paths.get(jobPaths.callInputsRoot.toString, callInputFile)

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
    val jobDescriptor: BackendJobDescriptor = jobDescriptorFromSingleCallWorkflow(workflowDescriptor, Map.empty, WorkflowOptions.empty, runtimeAttributeDefinitions)
    val backendRef = createBackendRef(jobDescriptor, emptyBackendConfig)
    val backend = backendRef.underlyingActor

    val execute = backend.execute
    backendRef ! AbortJobCommand

    whenReady(execute) { executionResponse =>
      executionResponse shouldBe a[AbortedResponse]
    }
  }

  def recoverSpec(completed: Boolean, writeReturnCode: Boolean = true): Assertion = {
    val workflowDescriptor = buildWorkflowDescriptor(HelloWorld)
    val jobDescriptor: BackendJobDescriptor = jobDescriptorFromSingleCallWorkflow(workflowDescriptor, Map.empty, WorkflowOptions.empty, runtimeAttributeDefinitions)
    val backendRef = createBackendRef(jobDescriptor, emptyBackendConfig)
    val backend = backendRef.underlyingActor

    val jobPaths = new JobPathsWithDocker(jobDescriptor.key, workflowDescriptor, ConfigFactory.empty)
    File(jobPaths.callExecutionRoot).createDirectories()
    File(jobPaths.stdout).write("Hello stubby ! ")
    File(jobPaths.stderr).touch()

    val pid =
      if (completed) {
        if (writeReturnCode)
          File(jobPaths.returnCode).write("0")
        "0"
      } else {
        import sys.process._
        val proc = Seq("bash", "-c", s"sleep 2; echo 0 > ${jobPaths.returnCode}").run()
        val pField = proc.getClass.getDeclaredField("p")
        pField.setAccessible(true)
        val p = pField.get(proc)
        val pidField = p.getClass.getDeclaredField("pid")
        pidField.setAccessible(true)
        pidField.get(p).toString
      }

    val execute = backend.recover

    val kvJobKey =
      KvJobKey(jobDescriptor.key.call.fullyQualifiedName, jobDescriptor.key.index, jobDescriptor.key.attempt)
    val scopedKey = ScopedKey(workflowDescriptor.id, kvJobKey, SharedFileSystemJob.JobIdKey)
    val kvPair = KvPair(scopedKey, Option(pid))

    backendRef ! kvPair

    whenReady(execute, Timeout(10.seconds.dilated)) { executionResponse =>
      if (writeReturnCode) {
        executionResponse should be(a[JobSucceededResponse])
        val succeededResponse = executionResponse.asInstanceOf[JobSucceededResponse]
        succeededResponse.returnCode.value should be(0)
        succeededResponse.jobOutputs should be(Map("salutation" -> JobOutput(WdlString("Hello stubby !"))))
      } else {
        executionResponse should be(a[JobFailedNonRetryableResponse])
        val failedResponse = executionResponse.asInstanceOf[JobFailedNonRetryableResponse]
        failedResponse.returnCode should be(empty)
        failedResponse.throwable should be(a[RuntimeException])
        failedResponse.throwable.getMessage should startWith("Unable to determine that 0 is alive, and")
        failedResponse.throwable.getMessage should endWith("call-hello/execution/rc does not exist.")
      }
    }
  }

  it should "recover a job in progress" in {
    recoverSpec(completed = false)
  }

  it should "recover a job that already completed" in {
    recoverSpec(completed = true)
  }

  it should "not recover a job for a non-existent pid" in {
    recoverSpec(completed = true, writeReturnCode = false)
  }

  it should "execute shards from a scatter" in {
    val workflowDescriptor = buildWorkflowDescriptor(TestWorkflows.Scatter)

    val call = workflowDescriptor.workflow.taskCalls.head

    0 to 2 foreach { shard =>
      // This assumes that engine will give us the evaluated value of the scatter item at the correct index
      // If this is not the case, more context/logic will need to be moved to the backend so it can figure it out by itself
      val symbolMaps: Map[LocallyQualifiedName, WdlInteger] = Map("scattering.intNumber" -> WdlInteger(shard))

      val runtimeAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, WorkflowOptions.empty)(call.task.runtimeAttributes.attrs)

      val jobDescriptor: BackendJobDescriptor =
        BackendJobDescriptor(workflowDescriptor, BackendJobDescriptorKey(call, Option(shard), 1), runtimeAttributes, fqnMapToDeclarationMap(symbolMaps))
      val backend = createBackend(jobDescriptor, emptyBackendConfig)
      val response =
        JobSucceededResponse(mock[BackendJobDescriptorKey], Some(0), Map("out" -> JobOutput(WdlInteger(shard))), None, Seq.empty)
      executeJobAndAssertOutputs(backend, response)
    }
  }

  it should "post process outputs" in {
    val inputFile = createCannedFile("localize", "content from json inputs").pathAsString
    val inputs = Map {
      "wf_localize.localize.inputFile" -> WdlFile(inputFile)
    }
    val workflowDescriptor = buildWorkflowDescriptor(OutputProcess, inputs)
    val jobDescriptor: BackendJobDescriptor = jobDescriptorFromSingleCallWorkflow(workflowDescriptor, inputs, WorkflowOptions.empty, runtimeAttributeDefinitions)
    val backend = createBackend(jobDescriptor, emptyBackendConfig)
    val jobPaths = new JobPathsWithDocker(jobDescriptor.key, workflowDescriptor, emptyBackendConfig.backendConfig)
    val expectedA = WdlFile(jobPaths.callExecutionRoot.resolve("a").toAbsolutePath.toString)
    val expectedB = WdlFile(jobPaths.callExecutionRoot.resolve("dir").toAbsolutePath.resolve("b").toString)
    val expectedOutputs = Map(
      "o1" -> JobOutput(expectedA),
      "o2" -> JobOutput(WdlArray(WdlArrayType(WdlFileType), Seq(expectedA, expectedB))),
      "o3" -> JobOutput(WdlFile(inputFile))
    )
    val expectedResponse = JobSucceededResponse(jobDescriptor.key, Some(0), expectedOutputs, None, Seq.empty)

    executeJobAndAssertOutputs(backend, expectedResponse)
  }

  it should "fail post processing if an output file is not found" in {
    val expectedResponse = JobFailedNonRetryableResponse(mock[BackendJobDescriptorKey],
      AggregatedException("Could not process output, file not found:", Seq.empty), Option(0))
    val workflow = TestWorkflow(buildWorkflowDescriptor(MissingOutputProcess), emptyBackendConfig, expectedResponse)
    val backend = createBackend(jobDescriptorFromSingleCallWorkflow(workflow.workflowDescriptor, Map.empty, WorkflowOptions.empty, runtimeAttributeDefinitions), workflow.config)
    testWorkflow(workflow, backend)
  }

  def createCannedFile(prefix: String, contents: String): File = {
    val suffix = ".out"
    File.newTemporaryFile(prefix, suffix).write(contents)
  }
}
