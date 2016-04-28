package cromwell.backend.impl.local

import java.nio.file.{Files, Paths}

import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionAbortedResponse, BackendJobExecutionFailedResponse, BackendJobExecutionSucceededResponse}
import cromwell.backend.BackendLifecycleActor.{BackendJobExecutionAbortFailedResponse, BackendJobExecutionAbortSucceededResponse}
import cromwell.backend.impl.local.BackendTestkitSpec.DockerTest
import cromwell.backend.impl.local.TestWorkflows._
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobDescriptorKey}
import cromwell.core._
import org.scalatest.FlatSpec
import org.scalatest.mock.MockitoSugar
import org.scalatest.prop.TableDrivenPropertyChecks
import wdl4s.WdlExpression
import wdl4s.expression.WdlFunctions
import wdl4s.types._
import wdl4s.values._

class LocalJobExecutionActorSpec extends FlatSpec with BackendTestkitSpec with MockitoSugar with TestFileUtil with TableDrivenPropertyChecks {

  behavior of "LocalBackend"

//  it should "execute an hello world workflow" in {
//    val expectedOutputs: CallOutputs = Map(
//      "salutation" -> CallOutput(WdlString("Hello you !"), None)
//    )
//    val expectedResponse = BackendJobExecutionSucceededResponse(mock[BackendJobDescriptorKey], expectedOutputs)
//    val wf = new TestWorkflow(buildWorkflowDescriptor(HelloWorld), defaultBackendConfig, expectedResponse)
//
//    testWorkflow(wf)
//  }
//
//  it should "execute an hello world workflow on Docker" taggedAs DockerTest in {
//    val expectedOutputs: CallOutputs = Map(
//      "salutation" -> CallOutput(WdlString("Hello you !"), None)
//    )
//    val expectedResponse = BackendJobExecutionSucceededResponse(mock[BackendJobDescriptorKey], expectedOutputs)
//    val wf = new TestWorkflow(buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" }"""), defaultBackendConfig, expectedResponse)
//
//    testWorkflow(wf)
//  }
//
//  it should "send back an execution failure if the task fails" in {
//    val expectedResponse = BackendJobExecutionFailedResponse(mock[BackendJobDescriptorKey], new Exception())
//    val wf = new TestWorkflow(buildWorkflowDescriptor(GoodbyeWorld), defaultBackendConfig, expectedResponse)
//
//    testWorkflow(wf)
//  }
//
//  it should "execute calls with input files and localize them appropriately" in {
//
//    def templateConf(localizers: String) = BackendConfigurationDescriptor("config",
//      ConfigFactory.parseString(
//        s"""
//           |config {
//           |  root = "local-cromwell-executions"
//           |  filesystems {
//           |    local {
//           |      localization = [
//           |        $localizers
//           |      ]
//           |    }
//           |  }
//           |}
//        """.stripMargin)
//    )
//
//    val hardConf = templateConf("hard-link")
//    val symConf = templateConf("soft-link")
//    val copyConf = templateConf("copy")
//
//    val jsonInputFile = createCannedFile("localize", "content from json inputs").toPath.toAbsolutePath.toString
//    val callInputFile = createCannedFile("localize", "content from call inputs").toPath.toAbsolutePath.toString
//    val inputs = Map(
//      "localize.workflowFile" -> WdlFile(callInputFile),
//      "localize.localize.inputFileFromJson" -> WdlFile(jsonInputFile)
//    )
//
//    val expectedOutputs: CallOutputs = Map(
//      "out" -> CallOutput(WdlArray(WdlArrayType(WdlStringType),
//        Array(
//          WdlString("content from json inputs"),
//          WdlString("content from call inputs"))
//      ), None)
//    )
//
//    val localizers = Table(
//      ("conf", "isSymLink"),
//      (hardConf, false),
//      (copyConf, false),
//      (symConf, true)
//    )
//
//    forAll(localizers) { (conf, isSymlink) =>
//      val wf = buildWorkflowDescriptor(InputFiles, inputs)
//      val backend = localBackend(jobDescriptorFromSingleCallWorkflow(wf), conf)
//      val jobDescriptor: BackendJobDescriptor = jobDescriptorFromSingleCallWorkflow(wf)
//      val expectedResponse = BackendJobExecutionSucceededResponse(jobDescriptor.key, expectedOutputs)
//
//      val jobPaths = new JobPaths(wf, defaultConfig, jobDescriptor.key)
//
//      whenReady(backend.execute) { executionResponse =>
//        assertResponse(executionResponse, expectedResponse)
//        val localizedJsonInputFile = Paths.get(jobPaths.callRoot.toString, jsonInputFile.toString)
//        val localizedCallInputFile = Paths.get(jobPaths.callRoot.toString, callInputFile.toString)
//
//        Files.isSymbolicLink(localizedJsonInputFile) shouldBe isSymlink
//        val realJsonInputFile = if (isSymlink) Files.readSymbolicLink(localizedJsonInputFile).toFile else localizedJsonInputFile.toFile
//        realJsonInputFile.exists() shouldBe true
//
//        Files.isSymbolicLink(localizedCallInputFile) shouldBe isSymlink
//        val realCallInputFile = if (isSymlink) Files.readSymbolicLink(localizedJsonInputFile).toFile else localizedCallInputFile.toFile
//        realCallInputFile.exists() shouldBe true
//      }
//    }
//  }

  it should "evaluate and coerce input values" in {
    val workflowDeclarations = Map(
      "str" -> WdlInteger(31380),
      "f1" -> WdlFloat(26)
    )

    val expectedOutputs: CallOutputs = Map(
      "outInts" -> CallOutput(WdlArray(WdlArrayType(WdlIntegerType), Array(WdlInteger(31380))), None),
      "outFloats" -> CallOutput(WdlArray(WdlArrayType(WdlFloatType), Array(WdlFloat(48), WdlFloat(63))), None)
    )
    val expectedResponse = BackendJobExecutionSucceededResponse(mock[BackendJobDescriptorKey], expectedOutputs)
    val wf = new TestWorkflow(buildWorkflowDescriptor(InputExpressions, workflowDeclarations), defaultBackendConfig, expectedResponse)

    testWorkflow(wf)
  }

//  it should "abort a job and kill a process" in {
//    val wf = buildWorkflowDescriptor(Sleep10)
//    val jobDescriptor: BackendJobDescriptor = jobDescriptorFromSingleCallWorkflow(wf)
//    val backend = localBackend(jobDescriptor, defaultBackendConfig)
//
//    val execute = backend.execute
//    val abort = backend.abortJob
//
//    whenReady(abort) { abortResponse =>
//      abortResponse shouldBe a[BackendJobExecutionAbortSucceededResponse]
//    }
//
//    whenReady(execute) { executionResponse =>
//      executionResponse shouldBe a[BackendJobExecutionAbortedResponse]
//    }
//  }
//
//  it should "fail to abort if execution hasn't been requested" in {
//    val wf = buildWorkflowDescriptor(Sleep10)
//    val jobDescriptor: BackendJobDescriptor = jobDescriptorFromSingleCallWorkflow(wf)
//    val backend = localBackend(jobDescriptor, defaultBackendConfig)
//
//    val abort = backend.abortJob
//
//    whenReady(abort) { abortResponse =>
//      abortResponse shouldBe a[BackendJobExecutionAbortFailedResponse]
//    }
//  }
//
//  it should "execute shards from a scatter" in {
//    val wf = buildWorkflowDescriptor(Scatter)
//
//    val call = wf.workflowNamespace.workflow.calls.head
//
//    0 to 2 foreach { shard =>
//      val evaluatorBuilder = buildEvaluatorBuilder(call, Map("i" -> WdlInteger(shard)))
//
//      val jd: BackendJobDescriptor = new BackendJobDescriptor(wf, new BackendJobDescriptorKey(call, Option(shard), 1), evaluatorBuilder, inputsFor(wf, call))
//      val backend = localBackend(jd, defaultBackendConfig)
//      val response = BackendJobExecutionSucceededResponse(mock[BackendJobDescriptorKey], Map("out" -> CallOutput(WdlInteger(shard), None)))
//      executeJobAndAssertOutputs(backend, response)
//    }
//  }
}
