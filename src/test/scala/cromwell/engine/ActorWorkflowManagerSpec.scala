package cromwell.engine

import akka.actor.ActorSystem
import akka.testkit.TestActorRef
import com.typesafe.config.ConfigFactory
import cromwell.CromwellSpec
import cromwell.HelloWorldActorSpec._
import cromwell.binding.FullyQualifiedName
import cromwell.binding.values.{WdlString, WdlValue}
import cromwell.engine.WorkflowManagerActor.{SubmitWorkflow, WorkflowOutputs, WorkflowStatus}
import cromwell.util.ActorUtil
import cromwell.util.SampleWdl.HelloWorld

import scala.language.{higherKinds, postfixOps, reflectiveCalls}


class ActorWorkflowManagerSpec extends CromwellSpec(ActorSystem("ActorWorkflowManagerSpec", ConfigFactory.parseString(Config))) {

  "An ActorWorkflowManager" should {
    "run the Hello World workflow" in {
      implicit val workflowManagerActor = TestActorRef(ActorWorkflowManager.props, self, "Test the ActorWorkflowManager")

      val workflowId = waitForHandledMessage(named = "Done") {
        ActorUtil.messageAndWait(SubmitWorkflow(HelloWorld.WdlSource, HelloWorld.RawInputs), _.mapTo[WorkflowId])
      }

      val status = ActorUtil.messageWaitAndGet(WorkflowStatus(workflowId), _.mapTo[Option[WorkflowState]])
      status shouldEqual WorkflowSucceeded

      val outputs = ActorUtil.messageWaitAndGet(WorkflowOutputs(workflowId), _.mapTo[Option[Map[FullyQualifiedName, WdlValue]]])

      val actual = outputs.map { case (k, WdlString(string)) => k -> string }
      actual shouldEqual Map(HelloWorld.OutputKey -> HelloWorld.OutputValue)
    }
  }
}
