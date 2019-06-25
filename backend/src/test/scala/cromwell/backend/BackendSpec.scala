package cromwell.backend

import _root_.wdl.draft2.model._
import _root_.wdl.transforms.draft2.wdlom2wom.WdlDraft2WomExecutableMakers._
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, JobFailedNonRetryableResponse, JobFailedRetryableResponse, JobSucceededResponse}
import cromwell.backend.io.TestWorkflows._
import cromwell.core.callcaching.NoDocker
import cromwell.core.labels.Labels
import cromwell.core.{HogGroup, WorkflowId, WorkflowOptions}
import common.exception.AggregatedException
import org.scalatest.Matchers
import org.scalatest.concurrent.{ScalaFutures, ScaledTimeSpans}
import org.scalatest.time.{Millis, Seconds, Span}
import org.specs2.mock.Mockito
import spray.json.{JsObject, JsValue}
import wom.callable.Callable.{InputDefinition, RequiredInputDefinition}
import wom.core.WorkflowSource
import wom.expression.{NoIoFunctionSet, WomExpression}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.{CommandCallNode, OptionalGraphInputNodeWithDefault}
import wom.values.WomValue
import wom.transforms.WomExecutableMaker.ops._

trait BackendSpec extends ScalaFutures with Matchers with Mockito with ScaledTimeSpans {

  implicit val defaultPatience = PatienceConfig(timeout = scaled(Span(10, Seconds)), interval = Span(500, Millis))

  def testWorkflow(workflow: TestWorkflow, backend: BackendJobExecutionActor, inputs: Map[String, WomValue] = Map.empty) = {
    executeJobAndAssertOutputs(backend, workflow.expectedResponse)
  }

  def buildWorkflowDescriptor(workflowSource: WorkflowSource,
                              inputFileAsJson: Option[String],
                              options: WorkflowOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue])),
                              runtime: String = "") = {
    val wdlNamespace = WdlNamespaceWithWorkflow.load(workflowSource.replaceAll("RUNTIME", runtime),
      Seq.empty[Draft2ImportResolver]).get
    val executable = wdlNamespace.toWomExecutable(inputFileAsJson, NoIoFunctionSet, strictValidation = true) match {
      case Left(errors) => fail(s"Fail to build wom executable: ${errors.toList.mkString(", ")}")
      case Right(e) => e
    }
    
    BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      executable.entryPoint,
      executable.resolvedExecutableInputs.flatMap({case (port, v) => v.select[WomValue] map { port -> _ }}),
      options,
      Labels.empty,
      HogGroup("foo"),
      List.empty,
      None
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
    val call = workflowDescriptor.callable.graph.nodes.collectFirst({ case t: CommandCallNode => t}).get
    val jobKey = BackendJobDescriptorKey(call, None, 1)
    
    val inputDeclarations: Map[InputDefinition, WomValue] = call.inputDefinitionMappings.map {
      case (inputDef, resolved) => inputDef ->
        resolved.select[WomValue].orElse(
          resolved.select[WomExpression]
            .map(
              _.evaluateValue(inputs, NoIoFunctionSet).getOrElse(fail("Can't evaluate input"))
            )
        ).orElse(
          resolved.select[OutputPort] flatMap {
            case known if workflowDescriptor.knownValues.contains(known) => Option(workflowDescriptor.knownValues(known))
            case hasDefault if hasDefault.graphNode.isInstanceOf[OptionalGraphInputNodeWithDefault] =>
              Option(hasDefault.graphNode.asInstanceOf[OptionalGraphInputNodeWithDefault].default
                .evaluateValue(inputs, NoIoFunctionSet).getOrElse(fail("Can't evaluate input")))
            case _ => None
          }
        ).getOrElse {
          inputs(inputDef.name)
        }
    }.toMap
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(call.callable.runtimeAttributes, NoIoFunctionSet, Map.empty).getOrElse(fail("Failed to evaluate runtime attributes")) // .get is OK here because this is a test
    val runtimeAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, options)(evaluatedAttributes)
    BackendJobDescriptor(workflowDescriptor, jobKey, runtimeAttributes, inputDeclarations, NoDocker, None, Map.empty)
  }

  def jobDescriptorFromSingleCallWorkflow(wdl: WorkflowSource,
                                          options: WorkflowOptions,
                                          runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition]): BackendJobDescriptor = {
    val workflowDescriptor = buildWdlWorkflowDescriptor(wdl)
    val call = workflowDescriptor.callable.graph.nodes.collectFirst({ case t: CommandCallNode => t}).get
    val jobKey = BackendJobDescriptorKey(call, None, 1)
    val inputDeclarations = fqnMapToDeclarationMap(workflowDescriptor.knownValues)
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(call.callable.runtimeAttributes, NoIoFunctionSet, inputDeclarations).getOrElse(fail("Failed to evaluate runtime attributes")) // .get is OK here because this is a test
    val runtimeAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, options)(evaluatedAttributes)
    BackendJobDescriptor(workflowDescriptor, jobKey, runtimeAttributes, inputDeclarations, NoDocker, None, Map.empty)
  }

  def jobDescriptorFromSingleCallWorkflow(wdl: WorkflowSource,
                                          runtime: String,
                                          attempt: Int,
                                          options: WorkflowOptions,
                                          runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition]): BackendJobDescriptor = {
    val workflowDescriptor = buildWdlWorkflowDescriptor(wdl, runtime = runtime)
    val call = workflowDescriptor.callable.graph.nodes.collectFirst({ case t: CommandCallNode => t}).get
    val jobKey = BackendJobDescriptorKey(call, None, attempt)
    val inputDeclarations = fqnMapToDeclarationMap(workflowDescriptor.knownValues)
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(call.callable.runtimeAttributes, NoIoFunctionSet, inputDeclarations).getOrElse(fail("Failed to evaluate runtime attributes")) // .get is OK here because this is a test
    val runtimeAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, options)(evaluatedAttributes)
    BackendJobDescriptor(workflowDescriptor, jobKey, runtimeAttributes, inputDeclarations, NoDocker, None, Map.empty)
  }

  def assertResponse(executionResponse: BackendJobExecutionResponse, expectedResponse: BackendJobExecutionResponse) = {
    (executionResponse, expectedResponse) match {
      case (JobSucceededResponse(_, _, responseOutputs, _, _, _, _), JobSucceededResponse(_, _, expectedOutputs, _, _, _, _)) =>
        responseOutputs.outputs.size shouldBe expectedOutputs.outputs.size
        responseOutputs.outputs foreach {
          case (fqn, out) =>
            val expectedOut = expectedOutputs.outputs.collectFirst({case (p, v) if p.name == fqn.name => v})
            expectedOut.getOrElse(fail(s"Output ${fqn.name} not found in ${expectedOutputs.outputs.map(_._1.name)}")).valueString shouldBe out.valueString
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
    val call = workflowDescriptor.callable.graph.nodes.collectFirst({ case t: CommandCallNode => t}).get
    BackendJobDescriptorKey(call, None, 1)
  }
}

object BackendSpec extends BackendSpec
