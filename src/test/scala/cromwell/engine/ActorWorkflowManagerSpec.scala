package cromwell.engine

import akka.actor.ActorSystem
import akka.testkit.TestActorRef
import com.typesafe.config.ConfigFactory
import cromwell.{binding, CromwellSpec}
import cromwell.HelloWorldActorSpec._
import cromwell.binding.FullyQualifiedName
import cromwell.binding.values.{WdlString, WdlValue}
import cromwell.engine.WorkflowManagerActor.{SubmitWorkflow, WorkflowOutputs, WorkflowStatus}
import cromwell.util.ActorTestUtil
import cromwell.util.SampleWdl.HelloWorld

import scala.language.{higherKinds, postfixOps, reflectiveCalls}


class ActorWorkflowManagerSpec extends CromwellSpec(ActorSystem("ActorWorkflowManagerSpec", ConfigFactory.parseString(Config))) {

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
