package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import cats.data.NonEmptyList
import com.dimafeng.testcontainers.Container
import cromwell.core.Tags.DbmsTest
import cromwell.core.{WorkflowId, WorkflowOptions, WorkflowSourceFilesCollection}
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.{WorkflowStoreAbortResponse, WorkflowStoreState}
import cromwell.services.database.{DatabaseSystem, DatabaseTestKit, EngineDatabaseType}
import org.scalatest.BeforeAndAfterAll
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.enablers.Emptiness._
import org.scalatest.time.{Millis, Seconds, Span}
import org.specs2.mock.Mockito

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

class SqlWorkflowStoreSpec extends AnyFlatSpec with Matchers with ScalaFutures with BeforeAndAfterAll with Mockito {
  implicit val ec = ExecutionContext.global
  implicit val defaultPatience = PatienceConfig(scaled(Span(10, Seconds)), scaled(Span(100, Millis)))
  val sourceFilesCollection = NonEmptyList.of(WorkflowSourceFilesCollection(Option("sample"), None, None, None, None, "input", WorkflowOptions.empty, "string", None, workflowOnHold = true, Seq.empty))

  DatabaseSystem.All foreach { databaseSystem =>

    behavior of s"SqlWorkflowStore on ${databaseSystem.name}"

    val containerOpt: Option[Container] = DatabaseTestKit.getDatabaseTestContainer(databaseSystem)

    lazy val dataAccess = DatabaseTestKit.initializeDatabaseByContainerOptTypeAndSystem(containerOpt, EngineDatabaseType, databaseSystem)

    lazy val workflowStore = SqlWorkflowStore(dataAccess)

    it should "start container if required" taggedAs DbmsTest in {
      containerOpt.foreach { _.start }
    }

    it should "honor the onHold flag" taggedAs DbmsTest in {
      (for {
        submissionResponses <- workflowStore.add(sourceFilesCollection)
        startableWorkflows <- workflowStore.fetchStartableWorkflows(10, "A00", 1.second)
        _ = startableWorkflows.map(_.id).intersect(submissionResponses.map(_.id).toList) should be(empty)
        _ <- workflowStore.switchOnHoldToSubmitted(submissionResponses.head.id)
        startableWorkflows2 <- workflowStore.fetchStartableWorkflows(10, "A00", 1.second)
        _ = startableWorkflows2.map(_.id).intersect(submissionResponses.map(_.id).toList).size should be(1)
      } yield ()).futureValue
    }

    it should "abort an onHold workflow" taggedAs DbmsTest in {
      (for {
        submissionResponses <- workflowStore.add(sourceFilesCollection)
        startableWorkflows <- workflowStore.fetchStartableWorkflows(10, "A01", 1.second)
        _ = startableWorkflows.map(_.id).intersect(submissionResponses.map(_.id).toList) should be(empty)
        abortWorkflowId = submissionResponses.head.id
        workflowStoreAbortResponse <- workflowStore.aborting(abortWorkflowId)
        _ = workflowStoreAbortResponse should be(WorkflowStoreAbortResponse.AbortedOnHoldOrSubmitted)
      } yield ()).futureValue
    }

    it should "abort an onHold then submitted workflow without a heartbeat" taggedAs DbmsTest in {
      (for {
        submissionResponses <- workflowStore.add(sourceFilesCollection)
        startableWorkflows <- workflowStore.fetchStartableWorkflows(10, "A02", 1.second)
        _ = startableWorkflows.map(_.id).intersect(submissionResponses.map(_.id).toList) should be(empty)
        abortWorkflowId = submissionResponses.head.id
        _ <- workflowStore.switchOnHoldToSubmitted(abortWorkflowId)
        workflowStoreAbortResponse <- workflowStore.aborting(abortWorkflowId)
        _ = workflowStoreAbortResponse should be(WorkflowStoreAbortResponse.AbortedOnHoldOrSubmitted)
      } yield ()).futureValue
    }

    it should "abort an onHold then submitted workflow with a heartbeat" taggedAs DbmsTest in {
      (for {
        submissionResponses <- workflowStore.add(sourceFilesCollection)
        startableWorkflows <- workflowStore.fetchStartableWorkflows(10, "A03", 1.second)
        _ = startableWorkflows.map(_.id).intersect(submissionResponses.map(_.id).toList) should be(empty)
        abortWorkflowId = submissionResponses.head.id
        _ <- workflowStore.switchOnHoldToSubmitted(abortWorkflowId)
        _ <- workflowStore.writeWorkflowHeartbeats(Set((abortWorkflowId, OffsetDateTime.now)), OffsetDateTime.now)
        workflowStoreAbortResponse <- workflowStore.aborting(abortWorkflowId)
        _ = workflowStoreAbortResponse should be(WorkflowStoreAbortResponse.AbortedOnHoldOrSubmitted)
      } yield ()).futureValue
    }

    it should "abort an onHold then running workflow without a heartbeat" taggedAs DbmsTest in {
      (for {
        submissionResponses <- workflowStore.add(sourceFilesCollection)
        startableWorkflows <- workflowStore.fetchStartableWorkflows(10, "A04", 1.second)
        _ = startableWorkflows.map(_.id).intersect(submissionResponses.map(_.id).toList) should be(empty)
        abortWorkflowId = submissionResponses.head.id
        _ <- workflowStore.switchOnHoldToSubmitted(abortWorkflowId)
        // Contact the db to change the status to running. Does not actually claim using the cromwellId.
        _ <- workflowStore.sqlDatabase.updateWorkflowState(
          abortWorkflowId.toString,
          WorkflowStoreState.Submitted.toString,
          WorkflowStoreState.Running.toString
        )
        workflowStoreAbortResponse <- workflowStore.aborting(abortWorkflowId)
        _ = workflowStoreAbortResponse should be(WorkflowStoreAbortResponse.AbortRequested)
      } yield ()).futureValue
    }

    it should "abort an onHold then running workflow with a heartbeat" taggedAs DbmsTest in {
      (for {
        submissionResponses <- workflowStore.add(sourceFilesCollection)
        startableWorkflows <- workflowStore.fetchStartableWorkflows(10, "A05", 1.second)
        _ = startableWorkflows.map(_.id).intersect(submissionResponses.map(_.id).toList) should be(empty)
        abortWorkflowId = submissionResponses.head.id
        _ <- workflowStore.switchOnHoldToSubmitted(abortWorkflowId)
        // Contact the db to change the status to running. Does not actually claim using the cromwellId.
        _ <- workflowStore.sqlDatabase.updateWorkflowState(
          abortWorkflowId.toString,
          WorkflowStoreState.Submitted.toString,
          WorkflowStoreState.Running.toString
        )
        _ <- workflowStore.writeWorkflowHeartbeats(Set((abortWorkflowId, OffsetDateTime.now)), OffsetDateTime.now)
        workflowStoreAbortResponse <- workflowStore.aborting(abortWorkflowId)
        _ = workflowStoreAbortResponse should be(WorkflowStoreAbortResponse.AbortRequested)
      } yield ()).futureValue
    }

    it should "not abort workflow that is not found" taggedAs DbmsTest in {
      val notFoundWorkflowId = WorkflowId.fromString("744e0645-1a1f-4ffe-a25d-a0be1f937fd7")

      (for {
        workflowStoreAbortResponse <- workflowStore.aborting(notFoundWorkflowId)
        _ = workflowStoreAbortResponse should be(WorkflowStoreAbortResponse.NotFound)
      } yield ()).futureValue
    }

    it should "stop container if required" taggedAs DbmsTest in {
      containerOpt.foreach { _.stop }
    }
  }
}
