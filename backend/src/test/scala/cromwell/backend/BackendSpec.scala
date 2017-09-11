package cromwell.backend

import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, JobFailedNonRetryableResponse, JobFailedRetryableResponse, JobSucceededResponse}
import cromwell.backend.io.TestWorkflows._
import cromwell.core.callcaching.NoDocker
import cromwell.core.labels.Labels
import cromwell.core.{NoIoFunctionSet, WorkflowId, WorkflowOptions}
import lenthall.exception.AggregatedException
import org.scalatest.Matchers
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.specs2.mock.Mockito
import spray.json.{JsObject, JsValue}
import wdl4s.wdl._
import wdl4s.wdl.values.WdlValue
import wdl4s.wom.callable.Callable.InputDefinition
import wdl4s.wom.graph.TaskCallNode

import scala.language.postfixOps

trait BackendSpec extends ScalaFutures with Matchers with Mockito {

  implicit val defaultPatience = PatienceConfig(timeout = Span(10, Seconds), interval = Span(500, Millis))

  def testWorkflow(workflow: TestWorkflow, backend: BackendJobExecutionActor, inputs: Map[String, WdlValue] = Map.empty) = {
    executeJobAndAssertOutputs(backend, workflow.expectedResponse)
  }

  def buildWorkflowDescriptor(workflowSource: WorkflowSource,
                              inputs: Map[String, WdlValue] = Map.empty,
                              options: WorkflowOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue])),
                              runtime: String = "") = {
    BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      WdlNamespaceWithWorkflow.load(workflowSource.replaceAll("RUNTIME", runtime),
        Seq.empty[ImportResolver]).get.workflow.womDefinition.getOrElse(fail("Cannot convert WdlWorkflow to WomDefinition")),
      inputs,
      options,
      Labels.empty
    )
  }

  def fqnMapToDeclarationMap(m: Map[String, WdlValue]): Map[InputDefinition, WdlValue] = {
    m map {
      case (fqn, v) =>
        val mockDeclaration = mock[InputDefinition]
        // TODO WOM: FIXME
//        mockDeclaration.fullyQualifiedName returns fqn
        mockDeclaration.name returns fqn.split('.').lastOption.getOrElse(fqn)
        mockDeclaration -> v
    }
  }

  def jobDescriptorFromSingleCallWorkflow(workflowDescriptor: BackendWorkflowDescriptor,
                                          inputs: Map[String, WdlValue],
                                          options: WorkflowOptions,
                                          runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition]): BackendJobDescriptor = {
    val call = workflowDescriptor.workflow.innerGraph.nodes.collectFirst({ case t: TaskCallNode => t}).get
    val jobKey = BackendJobDescriptorKey(call, None, 1)
    val inputDeclarations: Map[InputDefinition, WdlValue] = call.callable.inputs map { inputDefinition => inputDefinition -> inputs(inputDefinition.name) } toMap
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(call.callable.runtimeAttributes, NoIoFunctionSet, Map.empty).get // .get is OK here because this is a test
    val runtimeAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, options)(evaluatedAttributes)
    BackendJobDescriptor(workflowDescriptor, jobKey, runtimeAttributes, inputDeclarations, NoDocker, Map.empty)
  }

  def jobDescriptorFromSingleCallWorkflow(wdl: WorkflowSource,
                                          options: WorkflowOptions,
                                          runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition]): BackendJobDescriptor = {
    val workflowDescriptor = buildWorkflowDescriptor(wdl)
    val call = workflowDescriptor.workflow.innerGraph.nodes.collectFirst({ case t: TaskCallNode => t}).get
    val jobKey = BackendJobDescriptorKey(call, None, 1)
    val inputDeclarations = fqnMapToDeclarationMap(workflowDescriptor.knownValues)
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(call.callable.runtimeAttributes, NoIoFunctionSet, inputDeclarations).get // .get is OK here because this is a test
    val runtimeAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, options)(evaluatedAttributes)
    BackendJobDescriptor(workflowDescriptor, jobKey, runtimeAttributes, inputDeclarations, NoDocker, Map.empty)
  }

  def jobDescriptorFromSingleCallWorkflow(wdl: WorkflowSource,
                                          runtime: String,
                                          attempt: Int,
                                          options: WorkflowOptions,
                                          runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition]): BackendJobDescriptor = {
    val workflowDescriptor = buildWorkflowDescriptor(wdl, runtime = runtime)
    val call = workflowDescriptor.workflow.innerGraph.nodes.collectFirst({ case t: TaskCallNode => t}).get
    val jobKey = BackendJobDescriptorKey(call, None, attempt)
    val inputDeclarations = fqnMapToDeclarationMap(workflowDescriptor.knownValues)
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(call.callable.runtimeAttributes, NoIoFunctionSet, inputDeclarations).get // .get is OK here because this is a test
    val runtimeAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, options)(evaluatedAttributes)
    BackendJobDescriptor(workflowDescriptor, jobKey, runtimeAttributes, inputDeclarations, NoDocker, Map.empty)
  }

  def assertResponse(executionResponse: BackendJobExecutionResponse, expectedResponse: BackendJobExecutionResponse) = {
    (executionResponse, expectedResponse) match {
      case (JobSucceededResponse(_, _, responseOutputs, _, _, _), JobSucceededResponse(_, _, expectedOutputs, _, _, _)) =>
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

  def firstJobDescriptorKey(workflowDescriptor: BackendWorkflowDescriptor): BackendJobDescriptorKey = {
    val call = workflowDescriptor.workflow.innerGraph.nodes.collectFirst({ case t: TaskCallNode => t}).get
    BackendJobDescriptorKey(call, None, 1)
  }

  def firstJobDescriptor(workflowDescriptor: BackendWorkflowDescriptor,
                         inputs: Map[String, WdlValue] = Map.empty) = {
    BackendJobDescriptor(workflowDescriptor, firstJobDescriptorKey(workflowDescriptor), Map.empty, fqnMapToDeclarationMap(inputs), NoDocker, Map.empty)
  }
}

object BackendSpec extends BackendSpec
