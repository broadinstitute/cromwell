package cromwell.backend.io

import akka.actor.ActorSystem
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, FailedNonRetryableResponse, FailedRetryableResponse, SucceededResponse}
import cromwell.backend._
import cromwell.backend.io.TestWorkflows.TestWorkflow
import cromwell.core.{WorkflowId, WorkflowOptions}
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{Matchers, Tag}
import spray.json.{JsObject, JsValue}
import wdl4s._
import wdl4s.values.WdlValue

object BackendTestkitSpec {
  implicit val testActorSystem = ActorSystem("BackendTestActorSystem")
  object DockerTest extends Tag("DockerTest")
}

trait BackendTestkitSpec extends ScalaFutures with Matchers {

  implicit val defaultPatience = PatienceConfig(timeout = Span(5, Seconds), interval = Span(500, Millis))

  def testWorkflow(workflow: TestWorkflow, backend: BackendJobExecutionActor, inputs: Map[String, WdlValue] = Map.empty) = {
    executeJobAndAssertOutputs(backend, workflow.expectedResponse)
  }

  def buildWorkflowDescriptor(wdl: WdlSource,
                              inputs: Map[String, WdlValue] = Map.empty,
                              options: WorkflowOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue])),
                              runtime: String = "") = {
    new BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      NamespaceWithWorkflow.load(wdl.replaceAll("RUNTIME", runtime)),
      inputs,
      options
    )
  }

  def jobDescriptorFromSingleCallWorkflow(workflowDescriptor: BackendWorkflowDescriptor,
                                          inputs: Map[String, WdlValue] = Map.empty) = {
    val call = workflowDescriptor.workflowNamespace.workflow.calls.head
    val jobKey = new BackendJobDescriptorKey(call, None, 1)
    new BackendJobDescriptor(workflowDescriptor, jobKey, inputs)
  }

  def assertResponse(executionResponse: BackendJobExecutionResponse, expectedResponse: BackendJobExecutionResponse) = {
    (executionResponse, expectedResponse) match {
      case (SucceededResponse(_, responseOutputs), SucceededResponse(_, expectedOutputs)) =>
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

}
