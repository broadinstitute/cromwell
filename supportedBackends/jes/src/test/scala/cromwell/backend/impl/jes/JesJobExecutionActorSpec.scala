package cromwell.backend.impl.jes

import akka.actor.{Actor, ActorRef, Props}
import akka.testkit.{TestActorRef, TestProbe}
import cromwell.backend.BackendJobDescriptor
import cromwell.core.TestKitSuite
import org.scalatest.{FlatSpecLike, Matchers}
import org.specs2.mock.Mockito

import scala.concurrent.duration._
import akka.testkit._
import cromwell.backend.BackendJobExecutionActor.{ExecuteJobCommand, JobFailedNonRetryableResponse}
import cromwell.backend.impl.jes.ControllableFailingJabjea.JabjeaExplode

import scala.concurrent.{ExecutionContext, Promise}

class JesJobExecutionActorSpec extends TestKitSuite("JesJobExecutionActorSpec") with FlatSpecLike with Matchers with Mockito {

  behavior of "JesJobExecutionActor"

  private val AwaitAlmostNothing = 100.milliseconds.dilated
  private val TimeoutDuration = 10.seconds.dilated
  implicit val ec: ExecutionContext = system.dispatcher

  it should "catch failures in JABJEA initialization and fail the job accordingly" in {
    val jobDescriptor = mock[BackendJobDescriptor]
    val jesWorkflowInfo = mock[JesConfiguration]
    val initializationData = mock[JesBackendInitializationData]
    val serviceRegistryActor = system.actorOf(Props.empty)
    val jesBackendSingletonActor = Option(system.actorOf(Props.empty))

    val parent = TestProbe()
    val deathwatch = TestProbe()
    val params = JesSyncExecutionActorParams(jobDescriptor, jesWorkflowInfo, initializationData, serviceRegistryActor,
      jesBackendSingletonActor)
    val testJJEA = TestActorRef[TestJesJobExecutionActor](
      props = Props(new TestJesJobExecutionActor(params, Props(new ConstructorFailingJABJEA))),
      supervisor = parent.ref)
    deathwatch watch testJJEA

    // Nothing happens:
    parent.expectNoMsg(max = AwaitAlmostNothing)
    deathwatch.expectNoMsg(max = AwaitAlmostNothing)

    testJJEA.tell(msg = ExecuteJobCommand, sender = parent.ref)

    parent.expectMsgPF(max = TimeoutDuration) {
      case JobFailedNonRetryableResponse(_, throwable, _) =>
        throwable.getMessage should be("JesAsyncBackendJobExecutionActor failed and didn't catch its exception.")
    }
  }

  it should "catch failures at a random point during JABJEA processing and fail the job accordingly" in {
    val jobDescriptor = mock[BackendJobDescriptor]
    val jesWorkflowInfo = mock[JesConfiguration]
    val initializationData = mock[JesBackendInitializationData]
    val serviceRegistryActor = system.actorOf(Props.empty)
    val jesBackendSingletonActor = Option(system.actorOf(Props.empty))

    val parent = TestProbe()
    val deathwatch = TestProbe()
    val jabjeaConstructionPromise = Promise[ActorRef]()
    val params = JesSyncExecutionActorParams(jobDescriptor, jesWorkflowInfo, initializationData, serviceRegistryActor,
      jesBackendSingletonActor)
    val testJJEA = TestActorRef[TestJesJobExecutionActor](
      props = Props(new TestJesJobExecutionActor(params, Props(new ControllableFailingJabjea(jabjeaConstructionPromise)))),
      supervisor = parent.ref)
    deathwatch watch testJJEA

    // Nothing happens:
    parent.expectNoMsg(max = AwaitAlmostNothing)
    deathwatch.expectNoMsg(max = AwaitAlmostNothing)

    testJJEA.tell(msg = ExecuteJobCommand, sender = parent.ref)

    // Wait for the JABJEA to be spawned. Then kill it:
    parent.expectNoMsg(max = AwaitAlmostNothing)
    deathwatch.expectNoMsg(max = AwaitAlmostNothing)
    jabjeaConstructionPromise.future foreach { _ ! JabjeaExplode }

    parent.expectMsgPF(max = TimeoutDuration) {
      case JobFailedNonRetryableResponse(_, throwable, _) =>
        throwable.getMessage should be("JesAsyncBackendJobExecutionActor failed and didn't catch its exception.")
    }
  }
}

class TestJesJobExecutionActor(jesParams: JesSyncExecutionActorParams,
                               fakeJabjeaProps: Props) extends JesJobExecutionActor(jesParams) {
  override def jabjeaProps: Props = fakeJabjeaProps
}

class ConstructorFailingJABJEA extends ControllableFailingJabjea(Promise[ActorRef]()) {
  // Explode immediately in the constructor:
  explode()
}

class ControllableFailingJabjea(constructionPromise: Promise[ActorRef]) extends Actor {
  def explode(): Unit = {
    val boom = 1 == 1
    if (boom) throw new RuntimeException("Test Exception! Don't panic if this appears during a test run!")
  }
  constructionPromise.trySuccess(self)
  override def receive: Receive = {
    case JabjeaExplode => explode()
  }
}

object ControllableFailingJabjea {
  case object JabjeaExplode
}
