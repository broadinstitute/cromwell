package cromwell.backend.impl.jes

import akka.actor.{Actor, ActorRef, Props}
import akka.testkit._
import cromwell.backend.BackendJobExecutionActor.{ExecuteJobCommand, JobFailedNonRetryableResponse}
import cromwell.backend.impl.jes.ControllableFailingJabjea.JabjeaExplode
import cromwell.backend.standard.{DefaultStandardSyncExecutionActorParams, StandardSyncExecutionActor, StandardSyncExecutionActorParams}
import cromwell.backend.{BackendJobDescriptor, MinimumRuntimeSettings}
import cromwell.core.TestKitSuite
import org.scalatest.{FlatSpecLike, Matchers}
import org.specs2.mock.Mockito

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Promise}
import scala.util.{Failure, Success}

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
    val ioActor = system.actorOf(Props.empty)
    val jesBackendSingletonActor = Option(system.actorOf(Props.empty))

    initializationData.jesConfiguration returns jesWorkflowInfo

    val parent = TestProbe()
    val deathwatch = TestProbe()
    val params = DefaultStandardSyncExecutionActorParams(JesAsyncBackendJobExecutionActor.JesOperationIdKey, serviceRegistryActor, ioActor,
      jobDescriptor, null, Option(initializationData), jesBackendSingletonActor,
      classOf[JesAsyncBackendJobExecutionActor], MinimumRuntimeSettings())
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
    val ioActor = system.actorOf(Props.empty)
    val jesBackendSingletonActor = Option(system.actorOf(Props.empty))

    initializationData.jesConfiguration returns jesWorkflowInfo

    val parent = TestProbe()
    val deathwatch = TestProbe()
    val jabjeaConstructionPromise = Promise[ActorRef]()
    val params = DefaultStandardSyncExecutionActorParams(JesAsyncBackendJobExecutionActor.JesOperationIdKey, serviceRegistryActor, ioActor,
      jobDescriptor, null, Option(initializationData), jesBackendSingletonActor,
      classOf[JesAsyncBackendJobExecutionActor],
      MinimumRuntimeSettings())
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
    jabjeaConstructionPromise.future onComplete {
      case Success(jabjea) =>
        jabjea ! JabjeaExplode
      case Failure(throwable) =>
        val exception = new RuntimeException("Error creating jabjea for test!", throwable)
        exception.printStackTrace()
        throw exception
    }

    parent.expectMsgPF(max = TimeoutDuration) {
      case JobFailedNonRetryableResponse(_, throwable, _) =>
        throwable.getMessage should be("JesAsyncBackendJobExecutionActor failed and didn't catch its exception.")
    }
  }
}

class TestJesJobExecutionActor(params: StandardSyncExecutionActorParams,
                               fakeJabjeaProps: Props) extends StandardSyncExecutionActor(params) {
  override def createAsyncProps(): Props = fakeJabjeaProps
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
