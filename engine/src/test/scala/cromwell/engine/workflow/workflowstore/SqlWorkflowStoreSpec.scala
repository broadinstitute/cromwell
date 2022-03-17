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
import org.scalatest.time.{Millis, Seconds, Span}
import org.specs2.mock.Mockito
import spray.json.{JsObject, JsString}

import scala.concurrent.duration._
import scala.concurrent.{Await, ExecutionContext, Future}

class SqlWorkflowStoreSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with ScalaFutures with BeforeAndAfterAll with Mockito {
  implicit val ec = ExecutionContext.global
  implicit val defaultPatience = PatienceConfig(scaled(Span(20, Seconds)), scaled(Span(100, Millis)))

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

    def updateWfToRunning(startableWorkflows: List[WorkflowToStart]): Unit = {
      startableWorkflows.foreach { wf =>
        Await.result(workflowStore.sqlDatabase.updateWorkflowState(
          wf.id.toString,
          WorkflowStoreState.Submitted.toString,
          WorkflowStoreState.Running.toString
        ), 5.seconds)
      }
    }

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

    it should "start workflows from hog group with lowest count of running workflows" taggedAs DbmsTest in {
      // first submission of 50 workflows for hogGroup "Goldfinger"
      val goldFingerWorkflowIds = (for (_ <- 1 to 50) yield Await.result(workflowStore.add(includedGroupSourceFilesCollection1), 5.seconds)).flatMap(_.map(_.id).toList)

      // second submission of 50 workflows for hogGroup "Highlander"
      val highlanderWorkflowIds = (for (_ <- 1 to 50) yield Await.result(workflowStore.add(includedGroupSourceFilesCollection2), 5.seconds)).flatMap(_.map(_.id).toList)

      for (_ <- 1 to 10) yield {
        (for {
          // since both hog groups have 0 workflows running, the hog group with oldest submission time is picked first
          startableWorkflows1 <- workflowStore.fetchStartableWorkflows(5, "A08", 5.minutes, Set.empty[String])
          _ = startableWorkflows1.map(_.hogGroup.value).toSet.head should be("Goldfinger")
          _ = startableWorkflows1.map(_.id).foreach(x => goldFingerWorkflowIds.toList should contain(x))
          _ = updateWfToRunning(startableWorkflows1)

          startableWorkflows2 <- workflowStore.fetchStartableWorkflows(5, "A08", 5.minutes, Set.empty[String])
          _ = startableWorkflows2.map(_.hogGroup.value).toSet.head should be("Highlander")
          _ = startableWorkflows2.map(_.id).foreach(x => highlanderWorkflowIds.toList should contain(x))
          _ = updateWfToRunning(startableWorkflows2)
        } yield ()).futureValue
      }

      // remove entries from WorkflowStore
      (goldFingerWorkflowIds ++ highlanderWorkflowIds).foreach(id => Await.result(workflowStore.deleteFromStore(id), 5.seconds))
    }

    it should "respect excludedHogGroups and start workflows from hog group with lowest count of running workflows" taggedAs DbmsTest in {
      (for {
        // first submission of 10 workflows for hogGroup "Goldfinger"
        goldFingerSubmissions <- Future.sequence(for (_ <- 1 to 10) yield workflowStore.add(includedGroupSourceFilesCollection1))
        goldFingerWorkflowIds = goldFingerSubmissions.flatMap(_.map(_.id).toList)

        // second submission of 10 workflows for hogGroup "Zardoz"
        zardozSubmissions <- Future.sequence(for (_ <- 1 to 10) yield workflowStore.add(excludedGroupSourceFilesCollection))
        zardozWorkflowIds = zardozSubmissions.flatMap(_.map(_.id).toList)

        startableWorkflows1 <- workflowStore.fetchStartableWorkflows(5, "A08", 5.minutes, excludedGroups = Set("Zardoz"))
        _ = startableWorkflows1.map(_.hogGroup.value).toSet.head should be("Goldfinger")
        _ = startableWorkflows1.map(_.id).foreach(x => goldFingerWorkflowIds.toList should contain(x))
        _ = updateWfToRunning(startableWorkflows1)

        startableWorkflows2 <- workflowStore.fetchStartableWorkflows(5, "A08", 5.minutes, excludedGroups = Set("Zardoz"))
        _ = startableWorkflows2.map(_.hogGroup.value).toSet.head should be("Goldfinger")
        _ = startableWorkflows2.map(_.id).foreach(x => goldFingerWorkflowIds.toList should contain(x))
        _ = updateWfToRunning(startableWorkflows2)

        // there are 10 workflows from hog group "Zardoz" in the store, but since the group is excluded, 0 workflows are returned here
        startableWorkflows3 <- workflowStore.fetchStartableWorkflows(5, "A08", 5.minutes, excludedGroups = Set("Zardoz"))
        _ = startableWorkflows3.size should be(0)

        // hog group "Zardoz" has tokens to run workflows, hence don't exclude it
        startableWorkflows4 <- workflowStore.fetchStartableWorkflows(5, "A08", 5.minutes, Set.empty[String])
        _ = startableWorkflows4.map(_.hogGroup.value).toSet.head should be("Zardoz")
        _ = startableWorkflows4.map(_.id).foreach(x => zardozWorkflowIds.toList should contain(x))
        _ = updateWfToRunning(startableWorkflows4)

        startableWorkflows5 <- workflowStore.fetchStartableWorkflows(5, "A08", 5.minutes, Set.empty[String])
        _ = startableWorkflows5.map(_.hogGroup.value).toSet.head should be("Zardoz")
        _ = startableWorkflows5.map(_.id).foreach(x => zardozWorkflowIds.toList should contain(x))
        _ = updateWfToRunning(startableWorkflows5)

        // remove entries from WorkflowStore
        workflowsList = goldFingerWorkflowIds ++ zardozWorkflowIds
        _ = workflowsList.foreach(id => Await.result(workflowStore.deleteFromStore(id), 5.seconds))
      } yield()).futureValue
    }

    it should "start workflows from hog group with lowest count of running workflows for multiple hog groups" taggedAs DbmsTest in {
      (for {
        // first submission of 10 workflows for hogGroup "Goldfinger"
        goldFingerSubmissions <- Future.sequence(for (_ <- 1 to 10) yield workflowStore.add(includedGroupSourceFilesCollection1))
        goldFingerWorkflowIds = goldFingerSubmissions.flatMap(_.map(_.id).toList)

        // second submission of 10 workflows for hogGroup "Highlander"
        highlanderSubmissions <- Future.sequence(for (_ <- 1 to 15) yield workflowStore.add(includedGroupSourceFilesCollection2))
        highlanderWorkflowIds = highlanderSubmissions.flatMap(_.map(_.id).toList)

        // since both hog groups have 0 workflows running, the hog group with oldest submission time is picked first
        startableWorkflows1 <- workflowStore.fetchStartableWorkflows(5, "A08", 5.minutes, excludedGroups = Set.empty[String])
        _ = startableWorkflows1.map(_.hogGroup.value).toSet.head should be("Goldfinger")
        _ = startableWorkflows1.map(_.id).foreach(x => goldFingerWorkflowIds.toList should contain(x))
        _ = updateWfToRunning(startableWorkflows1)

        startableWorkflows2 <- workflowStore.fetchStartableWorkflows(5, "A08", 5.minutes, excludedGroups = Set.empty[String])
        _ = startableWorkflows2.map(_.hogGroup.value).toSet.head should be("Highlander")
        _ = startableWorkflows2.map(_.id).foreach(x => highlanderWorkflowIds.toList should contain(x))
        _ = updateWfToRunning(startableWorkflows2)

        // new submission for hog group "Finding Forrester"
        foresterSubmissions <- Future.sequence(for (_ <- 1 to 10) yield workflowStore.add(includedGroupSourceFilesCollection3))
        foresterWorkflowIds = foresterSubmissions.flatMap(_.map(_.id).toList)

        // now hog group "Finding Forrester" has 0 workflows running, hence it is picked to run
        startableWorkflows3 <- workflowStore.fetchStartableWorkflows(5, "A08", 5.minutes, excludedGroups = Set.empty[String])
        _ = startableWorkflows3.map(_.hogGroup.value).toSet.head should be("Finding Forrester")
        _ = startableWorkflows3.map(_.id).foreach(x => foresterWorkflowIds.toList should contain(x))
        _ = updateWfToRunning(startableWorkflows3)

        // since all 3 hog groups have 5 workflows running each, the hog group with oldest submission time is picked first
        startableWorkflows5 <- workflowStore.fetchStartableWorkflows(5, "A08", 5.minutes, excludedGroups = Set.empty[String])
        _ = startableWorkflows5.map(_.hogGroup.value).toSet.head should be("Goldfinger")
        _ = startableWorkflows5.map(_.id).foreach(x => goldFingerWorkflowIds.toList should contain(x))
        _ = updateWfToRunning(startableWorkflows5)

        // since both "Highlander" and "Finding Forrester" have 5 workflows in Running state, the hog group with oldest submission time is picked first
        startableWorkflows6 <- workflowStore.fetchStartableWorkflows(5, "A08", 5.minutes, excludedGroups = Set.empty[String])
        _ = startableWorkflows6.map(_.hogGroup.value).toSet.head should be("Highlander")
        _ = startableWorkflows6.map(_.id).foreach(x => highlanderWorkflowIds.toList should contain(x))
        _ = updateWfToRunning(startableWorkflows6)

        // "Finding Forrester" is now the hog group with least running workflows and has 5 more workflows to run, hence it is picked to run
        startableWorkflows4 <- workflowStore.fetchStartableWorkflows(5, "A08", 5.minutes, excludedGroups = Set.empty[String])
        _ = startableWorkflows4.map(_.hogGroup.value).toSet.head should be("Finding Forrester")
        _ = startableWorkflows4.map(_.id).foreach(x => foresterWorkflowIds.toList should contain(x))
        _ = updateWfToRunning(startableWorkflows4)

        startableWorkflows7 <- workflowStore.fetchStartableWorkflows(5, "A08", 5.minutes, excludedGroups = Set.empty[String])
        _ = startableWorkflows7.map(_.hogGroup.value).toSet.head should be("Highlander")
        _ = startableWorkflows7.map(_.id).foreach(x => highlanderWorkflowIds.toList should contain(x))
        _ = updateWfToRunning(startableWorkflows7)

        // remove entries from WorkflowStore
        workflowsList = goldFingerWorkflowIds ++ highlanderWorkflowIds ++ foresterWorkflowIds
        _ = workflowsList.foreach(id => Await.result(workflowStore.deleteFromStore(id), 5.seconds))
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
