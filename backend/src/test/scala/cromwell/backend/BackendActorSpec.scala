package cromwell.backend

import java.util.UUID

import akka.actor._
import akka.testkit.{DefaultTimeout, ImplicitSender, TestDuration, TestKit}
import cromwell.backend.BackendActor._
import cromwell.backend.TestWorkflowActor.{EventList, GetEvents, StopBackendActor}
import cromwell.backend.model.{JobDescriptor, WorkflowDescriptor}
import cromwell.core.eventbus.PubSubMediator
import org.scalatest.{BeforeAndAfter, BeforeAndAfterAll, Matchers, WordSpecLike}

import scala.collection.mutable.ArrayBuffer
import scala.concurrent.Future
import scala.concurrent.duration._

class BackendActorSpec extends TestKit(ActorSystem("BackendActorSpecSystem"))
  with DefaultTimeout with ImplicitSender with WordSpecLike with Matchers with PubSubMediator with BeforeAndAfter with BeforeAndAfterAll {

  val Timeout = 1.second.dilated

  var testWorkflowActor: ActorRef = _

  before {
    testWorkflowActor = system.actorOf(TestWorkflowActor.props())
  }

  after {
    system.stop(testWorkflowActor)
  }

  override def afterAll {
    TestKit.shutdownActorSystem(system)
  }

  "A Backend actor" should {
    "execute validate and beforeAll just after backend actor instance is created." in {
      within(Timeout) {
        testWorkflowActor ! GetEvents
        expectMsgPF() {
          case backendActorEvents: EventList =>
            if (!backendActorEvents.events.contains(ValidationSucceeded) || !backendActorEvents.events.contains(BeforeAllSucceeded)) {
              fail(s"All expected events were not raised. Expected events: ValidationSucceed and BeforeAllSucceeded.")
            }
          case unknown => fail(s"Response is not of type ArrayBuffer. Response: $unknown")
        }
      }
    }

    "execute afterAll just after backend actor instance is killed." in {
      within(Timeout) {
        testWorkflowActor ! StopBackendActor
        Thread.sleep(50) //Giving time to WorkflowActor to receive the event.
        testWorkflowActor ! GetEvents
        expectMsgPF() {
          case backendActorEvents: EventList =>
            if (!backendActorEvents.events.contains(ValidationSucceeded)
              || !backendActorEvents.events.contains(BeforeAllSucceeded)
              || !backendActorEvents.events.contains(AfterAllSucceeded)) {
              fail(s"All expected events were not raised. Expected events: ValidationSucceed, BeforeAllSucceeded and AfterAllSucceeded.")
            }
          case unknown => fail(s"Response is not of type ArrayBuffer. Response: $unknown")
        }
      }
    }
  }
}

object TestWorkflowActor {

  trait TestWorkflowActorMessage

  case object GetEvents extends TestWorkflowActorMessage

  case object StopBackendActor extends TestWorkflowActorMessage

  case class EventList(events: List[BackendActorEvent])

  def props(): Props = Props(new TestWorkflowActor(
    WorkflowDescriptor(UUID.fromString("8d772774-f76c-11e5-9ce9-5e5517507c66"))))

}

class TestWorkflowActor(workflowDescriptor: WorkflowDescriptor) extends Actor {
  val testBackendActor = context.actorOf(TestBackendActor.props(workflowDescriptor))
  var listOfEvents = ArrayBuffer[BackendActorEvent]()

  override def receive = {
    case GetEvents =>
      sender() ! EventList(listOfEvents.toList)
    case event: BackendActorEvent =>
      listOfEvents += event
    case StopBackendActor =>
      testBackendActor ! PoisonPill
  }
}

object TestBackendActor {

  def props(workflowDescriptor: WorkflowDescriptor): Props = Props(new TestBackendActor(
    workflowDescriptor))

}

class TestBackendActor(val workflowDescriptor: WorkflowDescriptor) extends BackendActor {
  /**
    * Tries to recover last status and information on an on going job in the backend.
    *
    * @return An RecoverEvent with the result.
    */
  override def recover(): RecoverEvent = ???

  /**
    * Executes a job base on a job description.
    *
    * @param jobDescriptor All information needed to execute a task in the backend.
    * @return An ExecutionEvent with the result.
    */
  override def execute(jobDescriptor: JobDescriptor): Future[ExecutionEvent] = ???

  /**
    * Stops all job execution.
    *
    * @return An AbortAllEvent with the result.
    */
  override def abortAll(): AbortAllEvent = ???

  /**
    * Executes validation on workflow descriptor in order to see if the workflow can be executed by the backend.
    */
  override def validate(): Unit = {
    context.parent ! ValidationSucceeded
  }

  /**
    * Stops a job execution.
    *
    * @param jobDescriptor All information related to a task.
    * @return An AbortEvent with the result.
    */
  override def abort(jobDescriptor: JobDescriptor): AbortEvent = ???

  /**
    * Registers code to be executed before the backend is ready for executing jobs for the specific workflow.
    */
  override def beforeAll(): Unit = {
    context.parent ! BeforeAllSucceeded
  }

  /**
    * Registers code to be executed after the backend finished executing all related tasks for the specific workflow.
    */
  override def afterAll(): Unit = {
    context.parent ! AfterAllSucceeded
  }
}
