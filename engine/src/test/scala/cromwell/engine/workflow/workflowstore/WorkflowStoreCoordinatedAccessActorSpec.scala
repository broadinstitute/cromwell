package cromwell.engine.workflow.workflowstore

import java.io.IOException
import java.time.OffsetDateTime

import akka.pattern._
import akka.testkit._
import akka.util.Timeout
import cats.data.NonEmptyVector
import cromwell.core._
import cromwell.engine.workflow.workflowstore.WorkflowStoreCoordinatedAccessActor.{FetchStartableWorkflows, WriteHeartbeats}
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{AsyncFlatSpecLike, Matchers}

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}

class WorkflowStoreCoordinatedAccessActorSpec extends TestKitSuite("WorkflowStoreCoordinatedWriteActorSpec")
  with AsyncFlatSpecLike with Matchers with TableDrivenPropertyChecks {

  behavior of "WorkflowStoreCoordinatedWriteActor"

  // So that we can timeout the asks below, change from the serial execution context to a parallel one
  override implicit def executionContext = scala.concurrent.ExecutionContext.Implicits.global

  def sleepAndThrow: Nothing = {
    Thread.sleep(30.seconds.dilated.toMillis)
    throw new RuntimeException("test should have timed out")
  }

  it should "writeHeartBeats" in {
    val expected = 12345
    val workflowStore = new InMemoryWorkflowStore {
      override def writeWorkflowHeartbeats(workflowIds: Set[(WorkflowId, OffsetDateTime)],
                                           heartbeatDateTime: OffsetDateTime)
                                          (implicit ec: ExecutionContext): Future[Int] = {
        Future.successful(expected)
      }
    }
    val actor = TestActorRef(new WorkflowStoreCoordinatedAccessActor(workflowStore))
    val request = WriteHeartbeats(NonEmptyVector.of((WorkflowId.randomId(), OffsetDateTime.now)), OffsetDateTime.now)
    implicit val timeout: Timeout = Timeout(2.seconds.dilated)
    actor.ask(request).mapTo[Int] map { actual =>
      actual should be(expected)
    }
  }

  it should "fetchStartableWorkflows" in {
    val collection = WorkflowSourceFilesCollection(
      workflowSource = Option("sample"),
      workflowUrl = None,
      workflowRoot = None,
      workflowType = None,
      workflowTypeVersion = None,
      inputsJson = "input",
      workflowOptions = WorkflowOptions.empty,
      labelsJson = "string",
      importsFile = None,
      workflowOnHold = true,
      warnings = Seq.empty
    )
    val now = OffsetDateTime.now()
    val expected: List[WorkflowToStart] = List(WorkflowToStart(WorkflowId.randomId(), now, collection, Submitted, HogGroup("foo")))
    val workflowStore = new InMemoryWorkflowStore {
      override def fetchStartableWorkflows(n: Int, cromwellId: String, heartbeatTtl: FiniteDuration)
                                          (implicit ec: ExecutionContext): Future[List[WorkflowToStart]] = {
        Future.successful(expected)
      }
    }
    val actor = TestActorRef(new WorkflowStoreCoordinatedAccessActor(workflowStore))
    val request = FetchStartableWorkflows(1, "test fetchStartableWorkflows success", 1.second)
    implicit val timeout: Timeout = Timeout(2.seconds.dilated)
    actor.ask(request).mapTo[List[WorkflowToStart]] map { actual =>
      actual should be(expected)
    }
  }

  it should "fetchStartableWorkflows with workflow url" in {
    val collection = WorkflowSourceFilesCollection(
      workflowSource = None,
      workflowUrl = Option("https://link-to-url"),
      workflowRoot = None,
      workflowType = None,
      workflowTypeVersion = None,
      inputsJson = "",
      workflowOptions = WorkflowOptions.empty,
      labelsJson = "",
      importsFile = None,
      workflowOnHold = false,
      warnings = Seq.empty
    )
    val now = OffsetDateTime.now()
    val expected: List[WorkflowToStart] = List(WorkflowToStart(WorkflowId.randomId(), now, collection, Submitted, HogGroup("foo")))
    val workflowStore = new InMemoryWorkflowStore {
      override def fetchStartableWorkflows(n: Int, cromwellId: String, heartbeatTtl: FiniteDuration)
                                          (implicit ec: ExecutionContext): Future[List[WorkflowToStart]] = {
        Future.successful(expected)
      }
    }
    val actor = TestActorRef(new WorkflowStoreCoordinatedAccessActor(workflowStore))
    val request = FetchStartableWorkflows(1, "test fetchStartableWorkflows with workflow url success", 1.second)
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
        override def writeWorkflowHeartbeats(workflowIds: Set[(WorkflowId, OffsetDateTime)],
                                             heartbeatDateTime: OffsetDateTime)
                                            (implicit ec: ExecutionContext): Future[Nothing] = {
          result()
        }
      }
      val actor = TestActorRef(new WorkflowStoreCoordinatedAccessActor(workflowStore))
      val request = WriteHeartbeats(NonEmptyVector.of((WorkflowId.randomId(), OffsetDateTime.now)), OffsetDateTime.now)
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
      val actor = TestActorRef(new WorkflowStoreCoordinatedAccessActor(workflowStore))
      val heartbeatTtlNotReallyUsed = 1.second
      val request = FetchStartableWorkflows(1, s"test $description fetchStartableWorkflows", heartbeatTtlNotReallyUsed)
      implicit val timeout: Timeout = Timeout(2.seconds.dilated)
      actor.ask(request).failed map { actual =>
        actual.getMessage should startWith(expectedMessagePrefix)
        actual.getClass should be(expectedException)
      }
    }
  }
}
