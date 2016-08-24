package cromwell.backend

import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, FailedNonRetryableResponse, FailedRetryableResponse, SucceededResponse}
import cromwell.backend.io.TestWorkflows._
import cromwell.core.{WorkflowId, WorkflowOptions}
import org.scalatest.Matchers
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import spray.json.{JsObject, JsValue}
import wdl4s._
import wdl4s.expression.NoFunctions
import wdl4s.values.WdlValue

trait BackendSpec extends ScalaFutures with Matchers {

  implicit val defaultPatience = PatienceConfig(timeout = Span(5, Seconds), interval = Span(500, Millis))

  def testWorkflow(workflow: TestWorkflow, backend: BackendJobExecutionActor, inputs: Map[String, WdlValue] = Map.empty) = {
    executeJobAndAssertOutputs(backend, workflow.expectedResponse)
  }

  def buildWorkflowDescriptor(wdl: WdlSource,
                              inputs: Map[String, WdlValue] = Map.empty,
                              options: WorkflowOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue])),
                              runtime: String = "") = {
    BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      NamespaceWithWorkflow.load(wdl.replaceAll("RUNTIME", runtime)),
      inputs,
      options
    )
  }

  def jobDescriptorFromSingleCallWorkflow(workflowDescriptor: BackendWorkflowDescriptor,
                                          inputs: Map[String, WdlValue],
                                          options: WorkflowOptions,
                                          runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition]): BackendJobDescriptor = {
    val call = workflowDescriptor.workflowNamespace.workflow.calls.head
    val jobKey = BackendJobDescriptorKey(call, None, 1)
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(call.task.runtimeAttributes, NoFunctions, inputs).get // .get is OK here because this is a test
    val runtimeAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, options)(evaluatedAttributes)
    BackendJobDescriptor(workflowDescriptor, jobKey, runtimeAttributes, inputs)
  }

  def jobDescriptorFromSingleCallWorkflow(wdl: WdlSource,
                                          options: WorkflowOptions,
                                          runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition]): BackendJobDescriptor = {
    val workflowDescriptor = buildWorkflowDescriptor(wdl)
    val call = workflowDescriptor.workflowNamespace.workflow.calls.head
    val jobKey = BackendJobDescriptorKey(call, None, 1)
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(call.task.runtimeAttributes, NoFunctions, workflowDescriptor.inputs).get // .get is OK here because this is a test
    val runtimeAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, options)(evaluatedAttributes)
    BackendJobDescriptor(workflowDescriptor, jobKey, runtimeAttributes, workflowDescriptor.inputs)
  }

  def jobDescriptorFromSingleCallWorkflow(wdl: WdlSource,
                                          runtime: String,
                                          attempt: Int,
                                          options: WorkflowOptions,
                                          runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition]): BackendJobDescriptor = {
    val workflowDescriptor = buildWorkflowDescriptor(wdl, runtime = runtime)
    val call = workflowDescriptor.workflowNamespace.workflow.calls.head
    val jobKey = BackendJobDescriptorKey(call, None, attempt)
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(call.task.runtimeAttributes, NoFunctions, workflowDescriptor.inputs).get // .get is OK here because this is a test
    val runtimeAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, options)(evaluatedAttributes)
    BackendJobDescriptor(workflowDescriptor, jobKey, runtimeAttributes, workflowDescriptor.inputs)
  }

  def assertResponse(executionResponse: BackendJobExecutionResponse, expectedResponse: BackendJobExecutionResponse) = {
    (executionResponse, expectedResponse) match {
      case (SucceededResponse(_, _, responseOutputs, _), SucceededResponse(_, _, expectedOutputs, _)) =>
        responseOutputs.size shouldBe expectedOutputs.size
        responseOutputs foreach {
          case (fqn, out) =>
            val expectedOut = expectedOutputs.get(fqn)
            expectedOut.isDefined shouldBe true
            expectedOut.get.wdlValue.valueString shouldBe out.wdlValue.valueString
        }
      case (FailedNonRetryableResponse(_, failure, _), FailedNonRetryableResponse(_, expectedFailure, _)) =>
        failure.getClass shouldBe expectedFailure.getClass
        failure.getMessage should include(expectedFailure.getMessage)
      case (FailedRetryableResponse(_, failure, _), FailedRetryableResponse(_, expectedFailure, _)) =>
        failure.getClass shouldBe expectedFailure.getClass
      case (response, expectation) =>
        fail(s"Execution response $response wasn't conform to expectation $expectation")
    }
  }

  def executeJobAndAssertOutputs(backend: BackendJobExecutionActor, expectedResponse: BackendJobExecutionResponse) = {
    whenReady(backend.execute) { executionResponse =>
      assertResponse(executionResponse, expectedResponse)
    }
  }

  lazy val emptyBackendConfig = BackendConfigurationDescriptor(
    ConfigFactory.parseString("{}"), ConfigFactory.load())

  def firstJobDescriptorKey(workflowDescriptor: BackendWorkflowDescriptor): BackendJobDescriptorKey = {
    val call = workflowDescriptor.workflowNamespace.workflow.calls.head
    BackendJobDescriptorKey(call, None, 1)
  }

  def firstJobDescriptor(workflowDescriptor: BackendWorkflowDescriptor,
                         inputs: Map[String, WdlValue] = Map.empty) = {
    BackendJobDescriptor(workflowDescriptor, firstJobDescriptorKey(workflowDescriptor), Map.empty, inputs)
  }
}

object BackendSpec extends BackendSpec
