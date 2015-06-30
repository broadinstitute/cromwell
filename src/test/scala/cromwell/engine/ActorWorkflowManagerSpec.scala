package cromwell.engine

import akka.testkit.TestActorRef
import cromwell.CromwellTestkitSpec
import cromwell.binding.values.WdlString
import cromwell.engine.WorkflowManagerActor.{SubmitWorkflow, WorkflowOutputs, WorkflowStatus}
import cromwell.util.ActorTestUtil
import cromwell.util.SampleWdl.HelloWorld
import cromwell.{CromwellSpec, binding}

import scala.language.{higherKinds, postfixOps, reflectiveCalls}

class ActorWorkflowManagerSpec extends CromwellTestkitSpec("ActorWorkflowManagerSpec") {
  "An ActorWorkflowManager" should {
    "run the Hello World workflow" in {
      implicit val workflowManagerActor = TestActorRef(WorkflowManagerActor.props, self, "Test the ActorWorkflowManager")

      val workflowId = waitForHandledMessagePattern(pattern = "Transition\\(.*,Running,Succeeded\\)$") {
        ActorTestUtil.messageAndWait(SubmitWorkflow(HelloWorld.WdlSource, HelloWorld.RawInputs), _.mapTo[WorkflowId])
      }

      val status = ActorTestUtil.messageAndWait(WorkflowStatus(workflowId), _.mapTo[Option[WorkflowState]]).get
      status shouldEqual WorkflowSucceeded

      val outputs = ActorTestUtil.messageAndWait(WorkflowOutputs(workflowId), _.mapTo[binding.WorkflowOutputs])
      val actual = outputs.map { case (k, WdlString(string)) => k -> string }
      actual shouldEqual Map(HelloWorld.OutputKey -> HelloWorld.OutputValue)
    }
  }
}
