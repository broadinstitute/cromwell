package cromwell.webservice

import akka.actor.{Actor, ActorRef, Props}
import akka.testkit.{TestActorRef, TestProbe}
import cromwell.core.TestKitSuite
import cromwell.webservice.EngineStatsActor.{EngineStats, JobCount, JobCountQuery}
import cromwell.webservice.EngineStatsActorSpec.FakeWorkflowActor
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._
import scala.language.postfixOps

class EngineStatsActorSpec extends TestKitSuite("EngineStatsActor") with FlatSpecLike with Matchers {
  behavior of "EngineStatsActor"

  val replyTo = TestProbe()
  val defaultTimeout = 500 millis

  it should "return double zeros with no WorkflowActors" in {
    TestActorRef(EngineStatsActor.props(List.empty[ActorRef], replyTo.ref, timeout = 200 millis))
    replyTo.expectMsg(defaultTimeout, EngineStats(0, 0))
  }

  it should "return snakeyes with a single workflow with one job" in {
    val workflowActors = List(Props(FakeWorkflowActor(1))) map { TestActorRef(_) }
    TestActorRef(EngineStatsActor.props(workflowActors, replyTo.ref, timeout = 200 millis))
    replyTo.expectMsg(defaultTimeout, EngineStats(1, 1))
  }

  it should "return an unemployed workflow when that's the world it lives in" in {
    val workflowActors = List(Props(FakeWorkflowActor(0))) map { TestActorRef(_) }
    TestActorRef(EngineStatsActor.props(workflowActors, replyTo.ref, timeout = 200 millis))
    replyTo.expectMsg(defaultTimeout, EngineStats(1, 0))
  }

  it should "return in a timely manner if one of the WorkflowActors is quiet" in {
    val workflowActors = List(Props(FakeWorkflowActor(0)), Props.empty) map { TestActorRef(_) }
    TestActorRef(EngineStatsActor.props(workflowActors, replyTo.ref, 50 millis))
    replyTo.expectMsg(defaultTimeout, EngineStats(1, 0))
  }

  it should "return the summation of jobs for all WorkflowActors" in {
    val workflowActors = List(Props(FakeWorkflowActor(1)), Props(FakeWorkflowActor(2))) map { TestActorRef(_) }
    TestActorRef(EngineStatsActor.props(workflowActors, replyTo.ref, timeout = 200 millis))
    replyTo.expectMsg(defaultTimeout, EngineStats(2, 3))
  }
}

object EngineStatsActorSpec {
  final case class FakeWorkflowActor(jobs: Int) extends Actor {
    override def receive = {
      case JobCountQuery => sender ! JobCount(jobs)
    }
  }
}
