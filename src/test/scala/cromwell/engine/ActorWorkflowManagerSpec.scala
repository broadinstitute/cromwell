package cromwell.engine

import akka.testkit.TestActorRef
import cromwell.binding.values.WdlString
import cromwell.engine.db.DummyDataAccess
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor.{SubmitWorkflow, WorkflowOutputs, WorkflowStatus}
import cromwell.util.SampleWdl.HelloWorld
import cromwell.{CromwellTestkitSpec, binding}

import scala.language.{higherKinds, postfixOps, reflectiveCalls}


class ActorWorkflowManagerSpec extends CromwellTestkitSpec("ActorWorkflowManagerSpec") {

  "An ActorWorkflowManager" should {
    "run the Hello World workflow" in {
      implicit val workflowManagerActor = TestActorRef(WorkflowManagerActor.props(DummyDataAccess()), self, "Test the ActorWorkflowManager")

      val workflowId = waitForHandledMessagePattern(pattern = "Transition\\(.*,Running,Succeeded\\)$") {
        messageAndWait[WorkflowId](SubmitWorkflow(HelloWorld.WdlSource, HelloWorld.RawInputs))
      }

      val status = messageAndWait[Option[WorkflowState]](WorkflowStatus(workflowId)).get
      status shouldEqual WorkflowSucceeded

      val outputs = messageAndWait[binding.WorkflowOutputs](WorkflowOutputs(workflowId))

      val actual = outputs.map { case (k, WdlString(string)) => k -> string }
      actual shouldEqual Map(HelloWorld.OutputKey -> HelloWorld.OutputValue)
    }
  }
}
