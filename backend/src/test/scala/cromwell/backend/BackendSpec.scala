package cromwell.backend

import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, JobFailedNonRetryableResponse, JobFailedRetryableResponse, JobSucceededResponse}
import cromwell.backend.io.TestWorkflows._
import cromwell.core.{WorkflowId, WorkflowOptions}
import lenthall.exception.AggregatedException
import org.scalatest.Matchers
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.specs2.mock.Mockito
import spray.json.{JsObject, JsValue}
import wdl4s._
import wdl4s.expression.NoFunctions
import wdl4s.values.WdlValue

trait BackendSpec extends ScalaFutures with Matchers with Mockito {

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
      WdlNamespaceWithWorkflow.load(wdl.replaceAll("RUNTIME", runtime), Seq.empty[ImportResolver]).workflow,
      inputs,
      options
    )
  }

  def fqnMapToDeclarationMap(m: Map[String, WdlValue]): Map[Declaration, WdlValue] = {
    m map {
      case (fqn, v) =>
        val mockDeclaration = mock[Declaration]
        mockDeclaration.fullyQualifiedName returns fqn
        mockDeclaration.unqualifiedName returns fqn.split('.').lastOption.getOrElse(fqn)
        mockDeclaration -> v
    }
  }

  def jobDescriptorFromSingleCallWorkflow(workflowDescriptor: BackendWorkflowDescriptor,
                                          inputs: Map[String, WdlValue],
                                          options: WorkflowOptions,
                                          runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition]): BackendJobDescriptor = {
    val call = workflowDescriptor.workflow.taskCalls.head
    val jobKey = BackendJobDescriptorKey(call, None, 1)
    val inputDeclarations = call.evaluateTaskInputs(inputs, NoFunctions)
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(call.task.runtimeAttributes, NoFunctions, inputDeclarations).get // .get is OK here because this is a test
    val runtimeAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, options)(evaluatedAttributes)
    BackendJobDescriptor(workflowDescriptor, jobKey, runtimeAttributes, inputDeclarations)
  }

  def jobDescriptorFromSingleCallWorkflow(wdl: WdlSource,
                                          options: WorkflowOptions,
                                          runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition]): BackendJobDescriptor = {
    val workflowDescriptor = buildWorkflowDescriptor(wdl)
    val call = workflowDescriptor.workflow.taskCalls.head
    val jobKey = BackendJobDescriptorKey(call, None, 1)
    val inputDeclarations = fqnMapToDeclarationMap(workflowDescriptor.knownValues)
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(call.task.runtimeAttributes, NoFunctions, inputDeclarations).get // .get is OK here because this is a test
    val runtimeAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, options)(evaluatedAttributes)
    BackendJobDescriptor(workflowDescriptor, jobKey, runtimeAttributes, inputDeclarations)
  }

  def jobDescriptorFromSingleCallWorkflow(wdl: WdlSource,
                                          runtime: String,
                                          attempt: Int,
                                          options: WorkflowOptions,
                                          runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition]): BackendJobDescriptor = {
    val workflowDescriptor = buildWorkflowDescriptor(wdl, runtime = runtime)
    val call = workflowDescriptor.workflow.taskCalls.head
    val jobKey = BackendJobDescriptorKey(call, None, attempt)
    val inputDeclarations = fqnMapToDeclarationMap(workflowDescriptor.knownValues)
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(call.task.runtimeAttributes, NoFunctions, inputDeclarations).get // .get is OK here because this is a test
    val runtimeAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, options)(evaluatedAttributes)
    BackendJobDescriptor(workflowDescriptor, jobKey, runtimeAttributes, inputDeclarations)
  }

  def assertResponse(executionResponse: BackendJobExecutionResponse, expectedResponse: BackendJobExecutionResponse) = {
    (executionResponse, expectedResponse) match {
      case (JobSucceededResponse(_, _, responseOutputs, _, _), JobSucceededResponse(_, _, expectedOutputs, _, _)) =>
        responseOutputs.size shouldBe expectedOutputs.size
        responseOutputs foreach {
          case (fqn, out) =>
            val expectedOut = expectedOutputs.get(fqn)
            expectedOut.isDefined shouldBe true
            expectedOut.get.wdlValue.valueString shouldBe out.wdlValue.valueString
        }
      case (JobFailedNonRetryableResponse(_, failure, _), JobFailedNonRetryableResponse(_, expectedFailure, _)) =>
        failure.getClass shouldBe expectedFailure.getClass
        concatenateCauseMessages(failure) should include(expectedFailure.getMessage)
      case (JobFailedRetryableResponse(_, failure, _), JobFailedRetryableResponse(_, expectedFailure, _)) =>
        failure.getClass shouldBe expectedFailure.getClass
      case (response, expectation) =>
        fail(s"Execution response $response wasn't conform to expectation $expectation")
    }
  }

  private def concatenateCauseMessages(t: Throwable): String = t match {
    case null => ""
    case ae: AggregatedException => ae.getMessage + " " + ae.throwables.map(innerT => concatenateCauseMessages(innerT.getCause)).mkString("\n")
    case other: Throwable => other.getMessage + concatenateCauseMessages(t.getCause)
  }


  def executeJobAndAssertOutputs(backend: BackendJobExecutionActor, expectedResponse: BackendJobExecutionResponse) = {
    whenReady(backend.execute) { executionResponse =>
      assertResponse(executionResponse, expectedResponse)
    }
  }

  lazy val emptyBackendConfig = BackendConfigurationDescriptor(
    ConfigFactory.parseString("{}"), ConfigFactory.load())

  def firstJobDescriptorKey(workflowDescriptor: BackendWorkflowDescriptor): BackendJobDescriptorKey = {
    val call = workflowDescriptor.workflow.taskCalls.head
    BackendJobDescriptorKey(call, None, 1)
  }

  def firstJobDescriptor(workflowDescriptor: BackendWorkflowDescriptor,
                         inputs: Map[String, WdlValue] = Map.empty) = {
    BackendJobDescriptor(workflowDescriptor, firstJobDescriptorKey(workflowDescriptor), Map.empty, fqnMapToDeclarationMap(inputs))
  }
}

object BackendSpec extends BackendSpec
