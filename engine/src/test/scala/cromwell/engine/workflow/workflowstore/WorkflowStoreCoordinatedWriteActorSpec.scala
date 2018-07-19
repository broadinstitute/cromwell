package cromwell.engine.workflow.workflowstore

import java.io.IOException

import akka.pattern._
import akka.testkit._
import akka.util.Timeout
import cats.data.NonEmptyVector
import cromwell.core.{TestKitSuite, WorkflowId, WorkflowSourceFilesCollection}
import cromwell.engine.workflow.workflowstore.WorkflowStoreCoordinatedWriteActor.{FetchStartableWorkflows, WriteHeartbeats}
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{AsyncFlatSpecLike, Matchers}

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}

class WorkflowStoreCoordinatedWriteActorSpec extends TestKitSuite("WorkflowStoreCoordinatedWriteActorSpec")
  with AsyncFlatSpecLike with Matchers with TableDrivenPropertyChecks {

  behavior of "WorkflowStoreCoordinatedWriteActor"

  // So that we can timeout the asks below, change from the serial execution context to a parallel one
  override implicit def executionContext = scala.concurrent.ExecutionContext.Implicits.global

  def sleepAndThrow: Nothing = {
    Thread.sleep(5 * 1000L)
    throw new RuntimeException("test should have timed out")
  }

  it should "writeHeartBeats" in {
    val expected = 12345
    val workflowStore = new InMemoryWorkflowStore {
      override def writeWorkflowHeartbeats(workflowIds: Set[WorkflowId])
                                          (implicit ec: ExecutionContext): Future[Int] = {
        Future.successful(expected)
      }
    }
    val actor = TestActorRef(new WorkflowStoreCoordinatedWriteActor(workflowStore))
    val request = WriteHeartbeats(NonEmptyVector.of(WorkflowId.randomId()))
    implicit val timeout: Timeout = Timeout(2.seconds.dilated)
    actor.ask(request).mapTo[Int] map { actual =>
      actual should be(expected)
    }
  }

  it should "fetchStartableWorkflows" in {
    val collection = WorkflowSourceFilesCollection(
      workflowSource = "sample",
      workflowUrl = None,
      workflowRoot = None,
      workflowType = None,
      workflowTypeVersion = None,
      inputsJson = "input",
      workflowOptionsJson = "option",
      labelsJson = "string",
      importsFile = None,
      workflowOnHold = true,
      warnings = Seq.empty
    )
    val expected: List[WorkflowToStart] = List(WorkflowToStart(WorkflowId.randomId(), collection, Submitted))
    val workflowStore = new InMemoryWorkflowStore {
      override def fetchStartableWorkflows(n: Int, cromwellId: String, heartbeatTtl: FiniteDuration)
                                          (implicit ec: ExecutionContext): Future[List[WorkflowToStart]] = {
        Future.successful(expected)
      }
    }
    val actor = TestActorRef(new WorkflowStoreCoordinatedWriteActor(workflowStore))
    val request = FetchStartableWorkflows(1, "test fetchStartableWorkflows success", 1.second)
    implicit val timeout: Timeout = Timeout(2.seconds.dilated)
    actor.ask(request).mapTo[List[WorkflowToStart]] map { actual =>
      actual should be(expected)
    }
  }

  val failureResponses = Table(
    ("description", "result", "expectedException", "expectedMessagePrefix"),
    ("a failure", () => Future.failed(new IOException("expected")), classOf[IOException], "expected"),
    ("a timeout", () => Future(sleepAndThrow), classOf[AskTimeoutException], "Ask timed out"),
  )

  forAll(failureResponses) { (description, result, expectedException, expectedMessagePrefix) =>
    it should s"fail to writeHeartBeats due to $description" in {
      val workflowStore = new InMemoryWorkflowStore {
        override def writeWorkflowHeartbeats(workflowIds: Set[WorkflowId])
                                            (implicit ec: ExecutionContext): Future[Nothing] = {
          result()
        }
      }
      val actor = TestActorRef(new WorkflowStoreCoordinatedWriteActor(workflowStore))
      val request = WriteHeartbeats(NonEmptyVector.of(WorkflowId.randomId()))
      implicit val timeout: Timeout = Timeout(2.seconds.dilated)
      actor.ask(request).failed map { actual =>
        actual.getMessage should startWith(expectedMessagePrefix)
        actual.getClass should be(expectedException)
      }
    }

    it should s"fail to fetchStartableWorkflows due to $description" in {
      val workflowStore = new InMemoryWorkflowStore {
        override def fetchStartableWorkflows(n: Int, cromwellId: String, heartbeatTtl: FiniteDuration)
                                            (implicit ec: ExecutionContext): Future[Nothing] = {
          result()
        }
      }
      val actor = TestActorRef(new WorkflowStoreCoordinatedWriteActor(workflowStore))
      val request = FetchStartableWorkflows(1, s"test $description fetchStartableWorkflows", 1.second)
      implicit val timeout: Timeout = Timeout(2.seconds.dilated)
      actor.ask(request).failed map { actual =>
        actual.getMessage should startWith(expectedMessagePrefix)
        actual.getClass should be(expectedException)
      }
    }
  }

}
