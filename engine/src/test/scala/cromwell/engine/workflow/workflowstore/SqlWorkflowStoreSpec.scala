package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import cats.data.NonEmptyList
import com.dimafeng.testcontainers.Container
import common.assertion.CromwellTimeoutSpec
import cromwell.core.Tags.DbmsTest
import cromwell.core.{WorkflowId, WorkflowOptions, WorkflowSourceFilesCollection, WorkflowSourceFilesWithoutImports}
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.{WorkflowStoreAbortResponse, WorkflowStoreState}
import cromwell.services.database.{DatabaseSystem, DatabaseTestKit, EngineDatabaseType, MetadataDatabaseType}
import org.scalatest.BeforeAndAfterAll
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.enablers.Emptiness._
import org.scalatest.time.{Millis, Seconds, Span}
import org.specs2.mock.Mockito
import spray.json.{JsObject, JsString}

import scala.concurrent.{ExecutionContext, Future}
import scala.concurrent.duration._

class SqlWorkflowStoreSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with ScalaFutures with BeforeAndAfterAll with Mockito {
  implicit val ec = ExecutionContext.global
  implicit val defaultPatience = PatienceConfig(scaled(Span(10, Seconds)), scaled(Span(100, Millis)))

  val onHoldSourceFilesCollection = NonEmptyList.of(
    WorkflowSourceFilesCollection(
      Option("sample"),
      None,
      None,
      None,
      None,
      "input",
      WorkflowOptions.empty,
      "string",
      None,
      workflowOnHold = true,
      Seq.empty,
      None
    )
  )

  val excludedGroupSourceFilesCollection = NonEmptyList.of(
    WorkflowSourceFilesCollection(
      Option("sample"),
      None,
      None,
      None,
      None,
      "input",
      WorkflowOptions(JsObject(Map("hogGroup" -> JsString("Zardoz")))),
      "string",
      None,
      workflowOnHold = false,
      Seq.empty,
      None
    )
  )

  val includedGroupSourceFilesCollection1 = NonEmptyList.of(
    WorkflowSourceFilesCollection(
      Option("sample"),
      None,
      None,
      None,
      None,
      "input",
      WorkflowOptions(JsObject(Map("hogGroup" -> JsString("Goldfinger")))),
      "string",
      None,
      workflowOnHold = false,
      Seq.empty,
      None
    )
  )

  val includedGroupSourceFilesCollection2 = NonEmptyList.of(
    WorkflowSourceFilesCollection(
      Option("sample"),
      None,
      None,
      None,
      None,
      "input",
      WorkflowOptions(JsObject(Map("hogGroup" -> JsString("Highlander")))),
      "string",
      None,
      workflowOnHold = false,
      Seq.empty,
      None
    )
  )

  val includedGroupSourceFilesCollection3 = NonEmptyList.of(
    WorkflowSourceFilesCollection(
      Option("sample"),
      None,
      None,
      None,
      None,
      "input",
      WorkflowOptions(JsObject(Map("hogGroup" -> JsString("Finding Forrester")))),
      "string",
      None,
      workflowOnHold = false,
      Seq.empty,
      None
    )
  )

  DatabaseSystem.All foreach { databaseSystem =>

    behavior of s"SqlWorkflowStore on ${databaseSystem.name}"

    val containerOpt: Option[Container] = DatabaseTestKit.getDatabaseTestContainer(databaseSystem)

    lazy val dataAccess = DatabaseTestKit.initializeDatabaseByContainerOptTypeAndSystem(containerOpt, EngineDatabaseType, databaseSystem)
    lazy val metadataDataAccess = DatabaseTestKit.initializeDatabaseByContainerOptTypeAndSystem(containerOpt, MetadataDatabaseType, databaseSystem)

    lazy val workflowStore = SqlWorkflowStore(dataAccess, metadataDataAccess)

    it should "start container if required" taggedAs DbmsTest in {
      containerOpt.foreach {
        _.start
      }
    }

    it should "honor the onHold flag" taggedAs DbmsTest in {
      (for {
        submissionResponses <- workflowStore.add(onHoldSourceFilesCollection)
        startableWorkflows <- workflowStore.fetchStartableWorkflows(10, "A00", 1.second, Set.empty)
        _ = startableWorkflows.map(_.id).intersect(submissionResponses.map(_.id).toList) should be(empty)
        _ <- workflowStore.switchOnHoldToSubmitted(submissionResponses.head.id)
        startableWorkflows2 <- workflowStore.fetchStartableWorkflows(10, "A00", 1.second, Set.empty)
        _ = startableWorkflows2.map(_.id).intersect(submissionResponses.map(_.id).toList).size should be(1)
        _ <- workflowStore.deleteFromStore(startableWorkflows2.head.id) // Tidy up
      } yield ()).futureValue
    }

    it should "abort an onHold workflow" taggedAs DbmsTest in {
      (for {
        submissionResponses <- workflowStore.add(onHoldSourceFilesCollection)
        startableWorkflows <- workflowStore.fetchStartableWorkflows(10, "A01", 1.second, Set.empty)
        _ = startableWorkflows.map(_.id).intersect(submissionResponses.map(_.id).toList) should be(empty)
        abortWorkflowId = submissionResponses.head.id
        workflowStoreAbortResponse <- workflowStore.abort(abortWorkflowId)
        _ = workflowStoreAbortResponse should be(WorkflowStoreAbortResponse.AbortedOnHoldOrSubmitted)
        _ <- workflowStore.deleteFromStore(abortWorkflowId) // Tidy up
      } yield ()).futureValue
    }

    it should "abort an onHold then submitted workflow without a heartbeat" taggedAs DbmsTest in {
      (for {
        submissionResponses <- workflowStore.add(onHoldSourceFilesCollection)
        startableWorkflows <- workflowStore.fetchStartableWorkflows(10, "A02", 1.second, Set.empty)
        _ = startableWorkflows.map(_.id).intersect(submissionResponses.map(_.id).toList) should be(empty)
        abortWorkflowId = submissionResponses.head.id
        _ <- workflowStore.switchOnHoldToSubmitted(abortWorkflowId)
        workflowStoreAbortResponse <- workflowStore.abort(abortWorkflowId)
        _ = workflowStoreAbortResponse should be(WorkflowStoreAbortResponse.AbortedOnHoldOrSubmitted)
        _ <- workflowStore.deleteFromStore(abortWorkflowId) // Tidy up
      } yield ()).futureValue
    }

    it should "abort an onHold then submitted workflow with a heartbeat" taggedAs DbmsTest in {
      (for {
        submissionResponses <- workflowStore.add(onHoldSourceFilesCollection)
        startableWorkflows <- workflowStore.fetchStartableWorkflows(10, "A03", 1.second, Set.empty)
        _ = startableWorkflows.map(_.id).intersect(submissionResponses.map(_.id).toList) should be(empty)
        abortWorkflowId = submissionResponses.head.id
        _ <- workflowStore.switchOnHoldToSubmitted(abortWorkflowId)
        _ <- workflowStore.writeWorkflowHeartbeats(Set((abortWorkflowId, OffsetDateTime.now)), OffsetDateTime.now)
        workflowStoreAbortResponse <- workflowStore.abort(abortWorkflowId)
        _ = workflowStoreAbortResponse should be(WorkflowStoreAbortResponse.AbortedOnHoldOrSubmitted)
        _ <- workflowStore.deleteFromStore(abortWorkflowId) // Tidy up
      } yield ()).futureValue
    }

    it should "abort an onHold then running workflow without a heartbeat" taggedAs DbmsTest in {
      (for {
        submissionResponses <- workflowStore.add(onHoldSourceFilesCollection)
        startableWorkflows <- workflowStore.fetchStartableWorkflows(10, "A04", 1.second, Set.empty)
        _ = startableWorkflows.map(_.id).intersect(submissionResponses.map(_.id).toList) should be(empty)
        abortWorkflowId = submissionResponses.head.id
        _ <- workflowStore.switchOnHoldToSubmitted(abortWorkflowId)
        // Contact the db to change the status to running. Does not actually claim using the cromwellId.
        _ <- workflowStore.sqlDatabase.updateWorkflowState(
          abortWorkflowId.toString,
          WorkflowStoreState.Submitted.toString,
          WorkflowStoreState.Running.toString
        )
        workflowStoreAbortResponse <- workflowStore.abort(abortWorkflowId)
        _ = workflowStoreAbortResponse should be(WorkflowStoreAbortResponse.AbortRequested)
        _ <- workflowStore.deleteFromStore(abortWorkflowId) // Tidy up
      } yield ()).futureValue
    }

    it should "abort an onHold then running workflow with a heartbeat" taggedAs DbmsTest in {
      (for {
        submissionResponses <- workflowStore.add(onHoldSourceFilesCollection)
        startableWorkflows <- workflowStore.fetchStartableWorkflows(10, "A05", 1.second, Set.empty)
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
        workflowStoreAbortResponse <- workflowStore.abort(abortWorkflowId)
        _ = workflowStoreAbortResponse should be(WorkflowStoreAbortResponse.AbortRequested)
        _ <- workflowStore.deleteFromStore(abortWorkflowId) // Tidy up
      } yield ()).futureValue
    }

    it should "not abort workflow that is not found" taggedAs DbmsTest in {
      val notFoundWorkflowId = WorkflowId.fromString("744e0645-1a1f-4ffe-a25d-a0be1f937fd7")

      (for {
        workflowStoreAbortResponse <- workflowStore.abort(notFoundWorkflowId)
        _ = workflowStoreAbortResponse should be(WorkflowStoreAbortResponse.NotFound)
      } yield ()).futureValue
    }

    it should "find a grouped workflow normally when not excluding" taggedAs DbmsTest in {
      (for {
        submissionResponses <- workflowStore.add(excludedGroupSourceFilesCollection)
        startableWorkflows <- workflowStore.fetchStartableWorkflows(10, "A06", 1.second, excludedGroups = Set.empty)
        _ = startableWorkflows.map(_.id).intersect(submissionResponses.map(_.id).toList).size should be(1)
        _ <- workflowStore.deleteFromStore(startableWorkflows.head.id) // Tidy up
      } yield ()).futureValue
    }

    it should "honor the excludedGroups parameter for a target group" taggedAs DbmsTest in {
      (for {
        submissionResponses <- workflowStore.add(excludedGroupSourceFilesCollection)
        startableWorkflows <- workflowStore.fetchStartableWorkflows(10, "A07", 1.second, excludedGroups = Set("Zardoz"))
        _ = startableWorkflows.map(_.id).intersect(submissionResponses.map(_.id).toList) should be(empty)
        _ <- workflowStore.deleteFromStore(submissionResponses.head.id) // Tidy up
      } yield ()).futureValue
    }

    it should "select appropriately with the excludedGroups parameter" taggedAs DbmsTest in {
      (for {
        submissionResponsesExcluded <- workflowStore.add(excludedGroupSourceFilesCollection)
        submissionResponsesIncluded1 <- workflowStore.add(includedGroupSourceFilesCollection1)
        submissionResponsesIncluded2 <- workflowStore.add(includedGroupSourceFilesCollection2)
        submissionResponsesIncluded3 <- workflowStore.add(includedGroupSourceFilesCollection3)
        startableWorkflows <- workflowStore.fetchStartableWorkflows(3, "A08", 1.second, excludedGroups = Set("Zardoz"))
        _ = startableWorkflows.map(_.id).intersect(submissionResponsesExcluded.map(_.id).toList).size should be(0)
        _ = startableWorkflows.map(_.id).intersect(submissionResponsesIncluded1.map(_.id).toList).size should be(1)
        _ = startableWorkflows.map(_.id).intersect(submissionResponsesIncluded2.map(_.id).toList).size should be(1)
        _ = startableWorkflows.map(_.id).intersect(submissionResponsesIncluded3.map(_.id).toList).size should be(1)
        _ = startableWorkflows.map(_.id).size should be(3)
        _ <- workflowStore.deleteFromStore(submissionResponsesExcluded.head.id) // Tidy up
        _ <- workflowStore.deleteFromStore(submissionResponsesIncluded1.head.id) // Tidy up
        _ <- workflowStore.deleteFromStore(submissionResponsesIncluded2.head.id) // Tidy up
        _ <- workflowStore.deleteFromStore(submissionResponsesIncluded3.head.id) // Tidy up
      } yield ()).futureValue
    }

    it should "accept and honor a requested workflow ID" taggedAs DbmsTest in {
      val requestedId = WorkflowId.randomId()

      val sourcesToSubmit = onHoldSourceFilesCollection.map(c => c.asInstanceOf[WorkflowSourceFilesWithoutImports].copy(
        requestedWorkflowId = Option(requestedId),
        workflowOnHold = false
      ))

      (for {
        submissionResponses <- workflowStore.add(sourcesToSubmit)
        startableWorkflows <- workflowStore.fetchStartableWorkflows(10, "A00", 1.second, Set.empty)
        _ = startableWorkflows.map(_.id).intersect(submissionResponses.map(_.id).toList).size should be(1)
        _ = startableWorkflows.map(_.id).intersect(submissionResponses.map(_.id).toList).head should be(requestedId)
        _ <- workflowStore.deleteFromStore(requestedId) // tidy up
      } yield ()).futureValue
    }

    it should "not accept a duplicate workflow ID" taggedAs DbmsTest in {
      val requestedId = WorkflowId.randomId()

      val workflowSourceFilesTemplate = onHoldSourceFilesCollection.head.asInstanceOf[WorkflowSourceFilesWithoutImports].copy(
        requestedWorkflowId = Option(requestedId)
      )

      val sourcesToSubmit1 = NonEmptyList.of(workflowSourceFilesTemplate)
      val sourcesToSubmit2 = NonEmptyList.of(workflowSourceFilesTemplate.copy(workflowOnHold = false))

      ((for {
        _ <- workflowStore.add(sourcesToSubmit1)
        _ <- workflowStore.add(sourcesToSubmit2)
      } yield ("incorrectly accepted")) recoverWith {
        case error => for {
          message <- Future {
            error.getMessage should be(s"Requested workflow IDs are already in use: $requestedId")
            "duplicate ID correctly detected"
          }
          stats <- workflowStore.stats
          _ = stats should be(Map(WorkflowStoreState.OnHold -> 1)) // Only the original (on-hold) version of requested ID 1 should be in the store
          _ <- workflowStore.deleteFromStore(requestedId) // tidy up
        } yield message
      }).futureValue should be("duplicate ID correctly detected")
    }

    it should "reject an entire workflow set if any requested workflow workflow ID is a duplicate" taggedAs DbmsTest in {
      val requestedId1 = WorkflowId.randomId()
      val requestedId2 = WorkflowId.randomId()
      val requestedId3 = WorkflowId.randomId()

      val workflowSourceFilesTemplate = onHoldSourceFilesCollection.head.asInstanceOf[WorkflowSourceFilesWithoutImports].copy(
        requestedWorkflowId = Option(requestedId1)
      )

      val sourcesToSubmit1 = NonEmptyList.of(workflowSourceFilesTemplate)

      val sourcesToSubmit2 = NonEmptyList.of(
        workflowSourceFilesTemplate.copy(requestedWorkflowId = Option(requestedId2), workflowOnHold = false),
        workflowSourceFilesTemplate.copy(requestedWorkflowId = Option(requestedId3), workflowOnHold = false),
        workflowSourceFilesTemplate.copy(requestedWorkflowId = Option(requestedId1), workflowOnHold = false) // duplicates the existing ID.
      )

      ((for {
        _ <- workflowStore.add(sourcesToSubmit1)
        _ <- workflowStore.add(sourcesToSubmit2)
      } yield ("incorrectly accepted")) recoverWith {
        case error => for {
          message <- Future {
            error.getMessage should be(s"Requested workflow IDs are already in use: $requestedId1")
            "duplicate ID correctly detected"
          }
          stats <- workflowStore.stats
          _ = stats should be(Map(WorkflowStoreState.OnHold -> 1)) // Only the original (on-hold) version of requested ID 1 should be in the store
          _ <- workflowStore.deleteFromStore(requestedId1)

        } yield message
      }).futureValue should be("duplicate ID correctly detected")
    }

    it should "reject an entire workflow set if it contains duplicate workflow ID within itself" taggedAs DbmsTest in {
      val requestedId1 = WorkflowId.randomId()
      val requestedId2 = WorkflowId.randomId()
      val requestedId3 = WorkflowId.randomId()

      val workflowSourceFilesTemplate = onHoldSourceFilesCollection.head.asInstanceOf[WorkflowSourceFilesWithoutImports]

      val sourcesToSubmit = NonEmptyList.of(
        workflowSourceFilesTemplate.copy(requestedWorkflowId = Option(requestedId1), workflowOnHold = false),
        workflowSourceFilesTemplate.copy(requestedWorkflowId = Option(requestedId2), workflowOnHold = false),
        workflowSourceFilesTemplate.copy(requestedWorkflowId = Option(requestedId3), workflowOnHold = false),
        workflowSourceFilesTemplate.copy(requestedWorkflowId = Option(requestedId1), workflowOnHold = false) // duplicates an ID already in the set
      )

      ((for {
        _ <- workflowStore.add(sourcesToSubmit)
      } yield ("incorrectly accepted")) recoverWith {
        case error => for {
          message <- Future {
            error.getMessage should be(s"Requested workflow IDs are duplicated: $requestedId1")
            "duplicate ID correctly detected"
          }
          stats <- workflowStore.stats
          _ = stats should be(Map.empty) // Nothing should have made it through to the store
          // No tidy up to do!
        } yield message
      }).futureValue should be("duplicate ID correctly detected")
    }

    it should "stop container if required" taggedAs DbmsTest in {
      containerOpt.foreach {
        _.stop
      }
    }
  }
}
