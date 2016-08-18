package cromwell.backend

import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, FailedNonRetryableResponse, FailedRetryableResponse, SucceededResponse}
import cromwell.backend.io.TestWorkflows._
import cromwell.core.{WorkflowId, WorkflowOptions}
import org.scalatest.Matchers
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import spray.json.{JsObject, JsValue}
import wdl4s.WdlExpression._
import wdl4s._
import wdl4s.expression.NoFunctions
import wdl4s.util.TryUtil
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
                                          inputs: Map[String, WdlValue] = Map.empty): BackendJobDescriptor = {
    val call = workflowDescriptor.workflowNamespace.workflow.calls.head
    val jobKey = BackendJobDescriptorKey(call, None, 1)
    BackendJobDescriptor(workflowDescriptor, jobKey, Map.empty, inputs)
  }

  def jobDescriptorFromSingleCallWorkflow(wdl: WdlSource): BackendJobDescriptor = {
    val workflowDescriptor = buildWorkflowDescriptor(wdl)
    val call = workflowDescriptor.workflowNamespace.workflow.calls.head
    val jobKey = BackendJobDescriptorKey(call, None, 1)
    BackendJobDescriptor(workflowDescriptor, jobKey, Map.empty, workflowDescriptor.inputs)
  }

  def jobDescriptorFromSingleCallWorkflow(wdl: WdlSource, runtime: String, attempt: Int): BackendJobDescriptor = {
    val workflowDescriptor = buildWorkflowDescriptor(wdl, runtime = runtime)
    val call = workflowDescriptor.workflowNamespace.workflow.calls.head
    val jobKey = BackendJobDescriptorKey(call, None, attempt)
    BackendJobDescriptor(workflowDescriptor, jobKey, Map.empty, workflowDescriptor.inputs)
  }

  def assertResponse(executionResponse: BackendJobExecutionResponse, expectedResponse: BackendJobExecutionResponse) = {
    (executionResponse, expectedResponse) match {
      case (SucceededResponse(_, _, responseOutputs), SucceededResponse(_, _, expectedOutputs)) =>
        responseOutputs.size shouldBe expectedOutputs.size
        responseOutputs foreach {
          case (fqn, out) =>
            val expectedOut = expectedOutputs.get(fqn)
            expectedOut.isDefined shouldBe true
            expectedOut.get.wdlValue.valueString shouldBe out.wdlValue.valueString
        }
      case (FailedNonRetryableResponse(_, failure, _), FailedNonRetryableResponse(_, expectedFailure, _)) =>
        // TODO improve this
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

  def createRuntimeAttributes(wdlSource: WdlSource, runtimeAttributes: String = "") = {
    val workflowDescriptor = buildWorkflowDescriptor(wdlSource, runtime = runtimeAttributes)

    def createLookup(call: Call): ScopedLookupFunction = {
      val declarations = workflowDescriptor.workflowNamespace.workflow.declarations ++ call.task.declarations
      val knownInputs = workflowDescriptor.inputs
      WdlExpression.standardLookupFunction(knownInputs, declarations, NoFunctions)
    }

    workflowDescriptor.workflowNamespace.workflow.calls map {
      call =>
        val ra = call.task.runtimeAttributes.attrs mapValues { _.evaluate(createLookup(call), NoFunctions) }
        TryUtil.sequenceMap(ra, "Runtime attributes evaluation").get
    }
  }
}

object BackendSpec extends BackendSpec
