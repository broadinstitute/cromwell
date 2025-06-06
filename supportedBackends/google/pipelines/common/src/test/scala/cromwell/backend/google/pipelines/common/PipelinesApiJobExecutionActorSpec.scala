package cromwell.backend.google.pipelines.common

import akka.actor.{Actor, ActorRef, Props}
import akka.testkit._
import cromwell.backend.BackendJobExecutionActor.{ExecuteJobCommand, JobFailedNonRetryableResponse}
import cromwell.backend.google.pipelines.common.ControllableFailingPabjea.JabjeaExplode
import cromwell.backend.standard.{
  DefaultStandardSyncExecutionActorParams,
  StandardSyncExecutionActor,
  StandardSyncExecutionActorParams
}
import cromwell.backend.{
  BackendJobDescriptor,
  BackendJobDescriptorKey,
  BackendWorkflowDescriptor,
  MinimumRuntimeSettings
}
import cromwell.core.{HogGroup, TestKitSuite, WorkflowId}
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import common.mock.MockSugar
import wom.graph.{CommandCallNode, WomIdentifier}

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Promise}
import scala.util.control.NoStackTrace
import scala.util.{Failure, Success}

class PipelinesApiJobExecutionActorSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with MockSugar {

  behavior of "PipelinesApiJobExecutionActor"

  private val AwaitAlmostNothing = 100.milliseconds.dilated
  private val TimeoutDuration = 10.seconds.dilated
  implicit val ec: ExecutionContext = system.dispatcher

  def jobDescriptor() =
    BackendJobDescriptor(
      BackendWorkflowDescriptor(WorkflowId.randomId(), null, Map.empty, null, null, HogGroup("asdf"), List.empty, None),
      BackendJobDescriptorKey(
        CommandCallNode(WomIdentifier.apply("asdf"), null, Set.empty, List.empty, Set.empty, null, None),
        None,
        0
      ),
      null,
      Map.empty,
      null,
      null,
      null
    )

  it should "catch failures in PABJEA initialization and fail the job accordingly" in {
    val jesWorkflowInfo = mock[PipelinesApiConfiguration]
    val initializationData = mock[PipelinesApiBackendInitializationData]
    val serviceRegistryActor = system.actorOf(Props.empty, "serviceRegistryActor-initialization")
    val ioActor = system.actorOf(Props.empty, "ioActor-initialization")
    val jesBackendSingletonActor = Option(system.actorOf(Props.empty, "jesBackendSingletonActor-initialization"))

    initializationData.papiConfiguration returns jesWorkflowInfo

    val parent = TestProbe("parent")
    val deathwatch = TestProbe("deathwatch")
    val params = DefaultStandardSyncExecutionActorParams(
      PipelinesApiAsyncBackendJobExecutionActor.JesOperationIdKey,
      serviceRegistryActor,
      ioActor,
      jobDescriptor(),
      null,
      Option(initializationData),
      jesBackendSingletonActor,
      classOf[PipelinesApiAsyncBackendJobExecutionActor],
      groupMetricsActor = system.actorOf(Props.empty, "groupMetricsActor-initialization"),
      MinimumRuntimeSettings()
    )
    val testJJEA = TestActorRef[TestPipelinesApiJobExecutionActor](
      props = Props(new TestPipelinesApiJobExecutionActor(params, Props(new ConstructorFailingJABJEA))),
      supervisor = parent.ref
    )
    deathwatch watch testJJEA

    // Nothing happens:
    parent.expectNoMessage(max = AwaitAlmostNothing)
    deathwatch.expectNoMessage(max = AwaitAlmostNothing)

    testJJEA.tell(msg = ExecuteJobCommand, sender = parent.ref)

    parent.expectMsgPF(max = TimeoutDuration) { case JobFailedNonRetryableResponse(_, throwable, _) =>
      throwable.getMessage should be(
        "PipelinesApiAsyncBackendJobExecutionActor failed and didn't catch its exception. This condition has been handled and the job will be marked as failed."
      )
    }
  }

  it should "catch failures at a random point during PABJEA processing and fail the job accordingly" in {
    val jesWorkflowInfo = mock[PipelinesApiConfiguration]
    val initializationData = mock[PipelinesApiBackendInitializationData]
    val serviceRegistryActor = system.actorOf(Props.empty, "serviceRegistryActor-random")
    val ioActor = system.actorOf(Props.empty, "ioActor-random")
    val jesBackendSingletonActor = Option(system.actorOf(Props.empty, "jesBackendSingletonActor-random"))

    initializationData.papiConfiguration returns jesWorkflowInfo

    val parent = TestProbe("parent")
    val deathwatch = TestProbe("deathwatch")
    val jabjeaConstructionPromise = Promise[ActorRef]()
    val params = DefaultStandardSyncExecutionActorParams(
      PipelinesApiAsyncBackendJobExecutionActor.JesOperationIdKey,
      serviceRegistryActor,
      ioActor,
      jobDescriptor(),
      null,
      Option(initializationData),
      jesBackendSingletonActor,
      classOf[PipelinesApiAsyncBackendJobExecutionActor],
      groupMetricsActor = system.actorOf(Props.empty, "groupMetricsActor-random"),
      MinimumRuntimeSettings()
    )
    val testJJEA = TestActorRef[TestPipelinesApiJobExecutionActor](
      props = Props(
        new TestPipelinesApiJobExecutionActor(params, Props(new ControllableFailingPabjea(jabjeaConstructionPromise)))
      ),
      supervisor = parent.ref
    )
    deathwatch watch testJJEA

    // Nothing happens:
    parent.expectNoMessage(max = AwaitAlmostNothing)
    deathwatch.expectNoMessage(max = AwaitAlmostNothing)

    testJJEA.tell(msg = ExecuteJobCommand, sender = parent.ref)

    // Wait for the JABJEA to be spawned. Then kill it:
    parent.expectNoMessage(max = AwaitAlmostNothing)
    deathwatch.expectNoMessage(max = AwaitAlmostNothing)
    jabjeaConstructionPromise.future onComplete {
      case Success(jabjea) =>
        jabjea ! JabjeaExplode
      case Failure(throwable) =>
        val exception = new RuntimeException("Error creating jabjea for test!", throwable)
        exception.printStackTrace()
        throw exception
    }

    parent.expectMsgPF(max = TimeoutDuration) { case JobFailedNonRetryableResponse(_, throwable, _) =>
      throwable.getMessage should be(
        "PipelinesApiAsyncBackendJobExecutionActor failed and didn't catch its exception. This condition has been handled and the job will be marked as failed."
      )
    }
  }
}

class TestPipelinesApiJobExecutionActor(params: StandardSyncExecutionActorParams, fakeJabjeaProps: Props)
    extends StandardSyncExecutionActor(params) {
  override def createAsyncProps(): Props = fakeJabjeaProps
}

class ConstructorFailingJABJEA extends ControllableFailingPabjea(Promise[ActorRef]()) {
  // Explode immediately in the constructor:
  explode()
}

class ControllableFailingPabjea(constructionPromise: Promise[ActorRef]) extends Actor {
  def explode(): Unit = {
    val boom = 1 == 1
    if (boom)
      throw new RuntimeException("Test Exception! Don't panic if this appears during a test run!") with NoStackTrace
  }
  constructionPromise.trySuccess(self)
  override def receive: Receive = { case JabjeaExplode =>
    explode()
  }
}

object ControllableFailingPabjea {
  case object JabjeaExplode
}
