package cromwell.engine

import akka.actor.{ActorRef, ActorSystem}
import akka.pattern.ask
import akka.testkit.TestActorRef
import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import cromwell.CromwellSpec
import cromwell.HelloWorldActorSpec._
import cromwell.binding.FullyQualifiedName
import cromwell.binding.values.{WdlString, WdlValue}
import cromwell.engine.WorkflowManagerActor.{SubmitWorkflow, WorkflowOutputs, WorkflowStatus}
import cromwell.util.SampleWdl.HelloWorld

import scala.concurrent.duration._
import scala.concurrent.{Await, Future}
import scala.language.{higherKinds, postfixOps, reflectiveCalls}


class ActorWorkflowManagerSpec extends CromwellSpec(ActorSystem("ActorWorkflowManagerSpec", ConfigFactory.parseString(Config))) {

  /**
   * Performs the following steps:
   *
   * <ol>
   * <li> Sends the specified message to the implicitly passed `ActorRef` via an `ask`.
   * <li> Collects the `Future[Any]` response.
   * <li> Invokes the downcasting function `downcast` on that `Future[Any]`.
   * <li> Issues a blocking `Await.result` on the `Future`.
   * <li> Invokes a `get` on the return of `Await.result`.
   * </ol>
   *
   * The return of `Await.result` is known to support a `get` method: the context bounds for `M` is the "duck type"
   * (formally, a structural type) having a `get` method returning a `U`.
   *
   */
  def messageWaitAndGet[U, M[U] <: {def get : U}](message: AnyRef, downcast: Future[_] => Future[M[U]])(implicit actorRef: ActorRef): U =
    messageAndWait(message, downcast)(actorRef).get

  /**
   * Performs the following steps:
   *
   * <ol>
   * <li> Sends the specified message to the implicitly passed `ActorRef` via an `ask`.
   * <li> Collects the `Future[Any]` response.
   * <li> Invokes the downcasting function `downcast` on that `Future[Any]`.
   * <li> Issues a blocking `Await.result` on the `Future`.
   * </ol>
   *
   */
  def messageAndWait[M](message: AnyRef, downcast: Future[_] => Future[M])(implicit actorRef: ActorRef): M = {
    val futureAny = actorRef ? message
    Await.result(downcast(futureAny), 5 seconds)
  }

  "An ActorWorkflowManager" should {
    "run the Hello World workflow" in {
      implicit val workflowManagerActor = TestActorRef(ActorWorkflowManager.props, self, "Test the ActorWorkflowManager")
      implicit val timeout = Timeout(5 seconds)

      val workflowId = waitForHandledMessage(named = "Done") {
        messageAndWait(SubmitWorkflow(HelloWorld.WdlSource, HelloWorld.RawInputs), _.mapTo[WorkflowId])
      }

      val status = messageWaitAndGet(WorkflowStatus(workflowId), _.mapTo[Option[WorkflowState]])
      status shouldEqual WorkflowSucceeded

      val outputs = messageWaitAndGet(WorkflowOutputs(workflowId), _.mapTo[Option[Map[FullyQualifiedName, WdlValue]]])

      val actual = outputs.map { case (k, WdlString(string)) => k -> string }
      actual shouldEqual Map(HelloWorld.OutputKey -> HelloWorld.OutputValue)
    }
  }
}
