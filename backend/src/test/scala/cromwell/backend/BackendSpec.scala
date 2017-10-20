package cromwell.backend

import _root_.wdl._
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, JobFailedNonRetryableResponse, JobFailedRetryableResponse, JobSucceededResponse}
import cromwell.backend.io.TestWorkflows._
import cromwell.core.callcaching.NoDocker
import cromwell.core.labels.Labels
import cromwell.core.{NoIoFunctionSet, WorkflowId, WorkflowOptions}
import common.exception.AggregatedException
import org.scalatest.Matchers
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.specs2.mock.Mockito
import spray.json.{JsObject, JsValue}
import wom.callable.Callable.{InputDefinition, RequiredInputDefinition}
import wom.core.WorkflowSource
import wom.expression.WomExpression
import wom.graph.GraphNodePort.OutputPort
import wom.graph.TaskCallNode
import wom.values.WomValue

trait BackendSpec extends ScalaFutures with Matchers with Mockito {

  implicit val defaultPatience = PatienceConfig(timeout = Span(10, Seconds), interval = Span(500, Millis))

  def testWorkflow(workflow: TestWorkflow, backend: BackendJobExecutionActor, inputs: Map[String, WomValue] = Map.empty) = {
    executeJobAndAssertOutputs(backend, workflow.expectedResponse)
  }

  def buildWorkflowDescriptor(workflowSource: WorkflowSource,
                              inputFileAsJson: Option[String],
                              options: WorkflowOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue])),
                              runtime: String = "") = {
    val wdlNamespace = WdlNamespaceWithWorkflow.load(workflowSource.replaceAll("RUNTIME", runtime),
      Seq.empty[ImportResolver]).get
    val executable = wdlNamespace.womExecutable(inputFileAsJson) match {
      case Left(errors) => fail(s"Fail to build wom executable: ${errors.toList.mkString(", ")}")
      case Right(e) => e
    }
    
    BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      wdlNamespace.workflow.womDefinition.getOrElse(fail("Cannot convert WdlWorkflow to WomDefinition")),
      executable.resolvedExecutableInputs.flatMap({case (port, v) => v.select[WomValue] map { port -> _ }}),
      options,
      Labels.empty
    )
  }

  def buildWdlWorkflowDescriptor(workflowSource: WorkflowSource,
                              inputFileAsJson: Option[String] = None,
                              options: WorkflowOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue])),
                              runtime: String = "") = {
    
    buildWorkflowDescriptor(workflowSource, inputFileAsJson, options, runtime)
  }

  def fqnWdlMapToDeclarationMap(m: Map[String, WomValue]): Map[InputDefinition, WomValue] = {
    m map {
      case (fqn, v) =>
        val mockDeclaration = RequiredInputDefinition(fqn, v.womType)
        mockDeclaration -> v
    }
  }

  def fqnMapToDeclarationMap(m: Map[OutputPort, WomValue]): Map[InputDefinition, WomValue] = {
    m map {
      case (outputPort, womValue) => RequiredInputDefinition(outputPort.name, womValue.womType) -> womValue 
    }
  }

  def jobDescriptorFromSingleCallWorkflow(workflowDescriptor: BackendWorkflowDescriptor,
                                          inputs: Map[String, WomValue],
                                          options: WorkflowOptions,
                                          runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition]): BackendJobDescriptor = {
    val call = workflowDescriptor.workflow.innerGraph.nodes.collectFirst({ case t: TaskCallNode => t}).get
    val jobKey = BackendJobDescriptorKey(call, None, 1)
    
    val inputDeclarations: Map[InputDefinition, WomValue] = call.inputDefinitionMappings.map {
      case (inputDef, resolved) => inputDef -> 
        resolved.select[WomValue].orElse(
          resolved.select[WomExpression]
            .map(
              _.evaluateValue(inputs, NoIoFunctionSet).getOrElse(fail("Can't evaluate input"))
            )
        ).orElse(
        workflowDescriptor.knownValues
          .get(resolved.select[OutputPort].get)
        )
        .getOrElse {
          inputs(inputDef.name) 
        }
    }
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(call.callable.runtimeAttributes, NoIoFunctionSet, Map.empty).getOrElse(fail("Failed to evaluate runtime attributes")) // .get is OK here because this is a test
    val runtimeAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, options)(evaluatedAttributes)
    BackendJobDescriptor(workflowDescriptor, jobKey, runtimeAttributes, inputDeclarations, NoDocker, Map.empty)
  }

  def jobDescriptorFromSingleCallWorkflow(wdl: WorkflowSource,
                                          options: WorkflowOptions,
                                          runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition]): BackendJobDescriptor = {
    val workflowDescriptor = buildWdlWorkflowDescriptor(wdl)
    val call = workflowDescriptor.workflow.innerGraph.nodes.collectFirst({ case t: TaskCallNode => t}).get
    val jobKey = BackendJobDescriptorKey(call, None, 1)
    val inputDeclarations = fqnMapToDeclarationMap(workflowDescriptor.knownValues)
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(call.callable.runtimeAttributes, NoIoFunctionSet, inputDeclarations).getOrElse(fail("Failed to evaluate runtime attributes")) // .get is OK here because this is a test
    val runtimeAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, options)(evaluatedAttributes)
    BackendJobDescriptor(workflowDescriptor, jobKey, runtimeAttributes, inputDeclarations, NoDocker, Map.empty)
  }

  def jobDescriptorFromSingleCallWorkflow(wdl: WorkflowSource,
                                          runtime: String,
                                          attempt: Int,
                                          options: WorkflowOptions,
                                          runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition]): BackendJobDescriptor = {
    val workflowDescriptor = buildWdlWorkflowDescriptor(wdl, runtime = runtime)
    val call = workflowDescriptor.workflow.innerGraph.nodes.collectFirst({ case t: TaskCallNode => t}).get
    val jobKey = BackendJobDescriptorKey(call, None, attempt)
    val inputDeclarations = fqnMapToDeclarationMap(workflowDescriptor.knownValues)
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(call.callable.runtimeAttributes, NoIoFunctionSet, inputDeclarations).getOrElse(fail("Failed to evaluate runtime attributes")) // .get is OK here because this is a test
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
            expectedOut.get.womValue.valueString shouldBe out.womValue.valueString
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
}

object BackendSpec extends BackendSpec
