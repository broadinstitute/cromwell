package cromwell.services.metadata.impl

import java.time.OffsetDateTime

import common.util.TimeUtil._
import cromwell.core.Tags.DbmsTest
import cromwell.core._
import cromwell.core.labels.{Label, Labels}
import cromwell.database.slick.MetadataSlickDatabase
import cromwell.services.MetadataServicesStore
import cromwell.services.database._
import cromwell.services.metadata._
import cromwell.services.metadata.impl.MetadataDatabaseAccess.SummaryResult
import org.scalatest.concurrent.PatienceConfiguration.{Interval, Timeout}
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import org.specs2.mock.Mockito

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}

object MetadataDatabaseAccessSpec {
  val AllowFalse = Seq(QueryParameter("allow", "false"))
  val AllowTrue = Seq(QueryParameter("allow", "true"))

  val Workflow1Name = "test1"
  val Workflow2Name = "test2"
  val Workflow3Name = "test3"

  val ParentWorkflowName = "test_parentWorkflow"
  val ParentWorkflow2Name = "test_parentWorkflow_2"
  val SubworkflowName = "test_subworkflow"
  val Subworkflow2Name = "test_subworkflow_2"
}

class MetadataDatabaseAccessSpec extends FlatSpec with Matchers with ScalaFutures with BeforeAndAfterAll with Mockito {
  import MetadataDatabaseAccessSpec._

  implicit val ec = ExecutionContext.global

  implicit val defaultPatience = PatienceConfig(scaled(Span(30, Seconds)), scaled(Span(100, Millis)))

  DatabaseSystem.All foreach { databaseSystem =>
    behavior of s"MetadataDatabaseAccess on ${databaseSystem.shortName}"

    lazy val dataAccess = new MetadataDatabaseAccess with MetadataServicesStore {
      override val metadataDatabaseInterface: MetadataSlickDatabase = {
        // NOTE: EngineLiquibaseSettings **MUST** always run before the MetadataLiquibaseSettings
        DatabaseTestKit.initializedDatabaseFromSystem(EngineDatabaseType, databaseSystem)
        DatabaseTestKit.initializedDatabaseFromSystem(MetadataDatabaseType, databaseSystem)
      }
    }

    def publishMetadataEvents(baseKey: MetadataKey, keyValues: Array[(String, String)]): Future[Unit] = {
      val events = keyValues map { case (k, v) =>
        MetadataEvent(baseKey.copy(key = k), MetadataValue(v))
      }
      dataAccess.addMetadataEvents(events)
    }

    def baseWorkflowMetadata(name: String, labels: Set[Label] = Set.empty): Future[WorkflowId] = {
      val workflowId = WorkflowId.randomId()
      val defaultLabels = Set(Label("cromwell-workflow-name", name))
      val labelMetadata = (labels ++ defaultLabels).map(label => (s"${WorkflowMetadataKeys.Labels}:${label.key}", label.value)).toArray

      val workflowKey = MetadataKey(workflowId, jobKey = None, key = null)
      def keyAndValue(name: String) = Array(
        (WorkflowMetadataKeys.SubmissionTime, OffsetDateTime.now.toUtcMilliString),
        (WorkflowMetadataKeys.Status, WorkflowSubmitted.toString),
        (WorkflowMetadataKeys.Name, name),
        (WorkflowMetadataKeys.StartTime, OffsetDateTime.now.toUtcMilliString)
      ) ++ labelMetadata

      publishMetadataEvents(workflowKey, keyAndValue(name)).map(_ => workflowId)
    }

    def subworkflowMetadata(parentWorkflowId: WorkflowId, subworkflowName: String): Future[WorkflowId] = {
      val workflowId = WorkflowId.randomId()
      val workflowKey = MetadataKey(workflowId, jobKey = None, key = null)
      val metadataKeys = Array(
        (WorkflowMetadataKeys.Status, WorkflowRunning.toString),
        (WorkflowMetadataKeys.Name, subworkflowName),
        (WorkflowMetadataKeys.StartTime, OffsetDateTime.now.toUtcMilliString),
        (WorkflowMetadataKeys.ParentWorkflowId, parentWorkflowId.toString),
        (WorkflowMetadataKeys.RootWorkflowId, parentWorkflowId.toString),
      )

      publishMetadataEvents(workflowKey, metadataKeys).map(_ => workflowId)
    }

    it should "return pagination metadata only when page and pagesize query params are specified" taggedAs DbmsTest in {
      (for {
        _ <- baseWorkflowMetadata(Workflow1Name)
        //get metadata when page and pagesize are specified
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.Page.name -> "1", WorkflowQueryKey.PageSize.name -> "50"))) map { case (_, meta) =>
          meta match {
            case Some(_) =>
            case None => fail("Should have metadata when page and pagesize are specified.")
          }
        }
        //don't get metadata when page and pagesize are not specified
        _ <- dataAccess.queryWorkflowSummaries(
          WorkflowQueryParameters(Seq())) map { case (_, meta) =>
          meta match {
            case Some(_) => fail("Should not have metadata when page and pagesize are not specified")
            case None =>
          }
        }
      } yield()).futureValue
    }

    it should "sort metadata events by timestamp from older to newer" taggedAs DbmsTest in {

      def unorderedEvents(id: WorkflowId): Future[Vector[MetadataEvent]] = {
        val workflowKey = MetadataKey(id, jobKey = None, key = null)
        val now = OffsetDateTime.now()
        val yesterday = now.minusDays(1)
        val tomorrow = now.plusDays(1)

        val yesterdayEvent = MetadataEvent(workflowKey.copy(key = WorkflowMetadataKeys.WorkflowRoot), Option(MetadataValue("A")), yesterday)
        val nowEvent = MetadataEvent(workflowKey.copy(key = WorkflowMetadataKeys.WorkflowRoot), Option(MetadataValue("B")), now)
        val tomorrowEvent = MetadataEvent(workflowKey.copy(key = WorkflowMetadataKeys.WorkflowRoot), Option(MetadataValue("C")), tomorrow)

        val events = Vector(tomorrowEvent, yesterdayEvent, nowEvent)

        val expectedEvents = Vector(yesterdayEvent, nowEvent, tomorrowEvent)

        dataAccess.addMetadataEvents(events) map { _ => expectedEvents }
      }

      (for {
        workflow1Id <- baseWorkflowMetadata(Workflow1Name)
        expected <- unorderedEvents(workflow1Id)
        response <- dataAccess.queryMetadataEvents(MetadataQuery(workflow1Id, None, Option(WorkflowMetadataKeys.WorkflowRoot), None, None, expandSubWorkflows = false), 5.seconds)
        _ = response shouldBe expected
      } yield()).futureValue
    }

    def assertRowsProcessedAndSummarizationComplete(summaryResult: SummaryResult) = {
      withClue(s"asserting correctness of $summaryResult") {
        summaryResult.rowsProcessedIncreasing should be > 0L
        summaryResult.rowsProcessedDecreasing should be(0L)

        summaryResult.increasingGap should be(0L)
        summaryResult.decreasingGap should be(0L)
      }
    }

    it should "create and query a workflow" taggedAs DbmsTest in {

      val randomIds = Seq.fill(10)(WorkflowId.randomId().toString)

      val testLabel1 = Label("testing-key-1", "testing-value-1")
      val testLabel2 = Label("testing-key-2", "testing-value-2")
      val testLabel3 = Label("testing-key-3", "testing-value-3")

      def succeededWorkflowMetadata(id: WorkflowId): Future[Unit] = {
        val workflowKey = MetadataKey(id, jobKey = None, key = null)
        val keyAndValue = Array(
          (WorkflowMetadataKeys.Status, WorkflowRunning.toString),
          (WorkflowMetadataKeys.Status, WorkflowSucceeded.toString),
          (WorkflowMetadataKeys.EndTime, OffsetDateTime.now.toUtcMilliString))

        publishMetadataEvents(workflowKey, keyAndValue)
      }

      (for {
        workflow1Id <- baseWorkflowMetadata(Workflow1Name, Set(testLabel1, testLabel2))
        _ <- succeededWorkflowMetadata(workflow1Id)
        // Put a bit of space between the two workflows
        _ = Thread.sleep(50)
        workflow2Id <- baseWorkflowMetadata(Workflow2Name, Set(testLabel2, testLabel3))

        // refresh the metadata
        _ <- dataAccess.refreshWorkflowMetadataSummaries(1000) map assertRowsProcessedAndSummarizationComplete

        // Query with no filters
        (workflowQueryResult, workflowQueryResult2) <-
        dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq.empty)) map { case (response, _) =>
          val result = response.results find { r => r.name.contains(Workflow1Name) && r.end.isDefined } getOrElse
            fail(s"$Workflow1Name with an end not found in ${response.results}")
          val result2 = response.results find {
            _.name.contains(Workflow2Name)
          } getOrElse fail(s"$Workflow2Name not found in ${response.results}")
          (result, result2)
        }
        // Filter by name
        _ <- dataAccess.queryWorkflowSummaries(
          WorkflowQueryParameters(Seq(WorkflowQueryKey.Name.name -> Workflow1Name))
        ) map { case (response, _) =>
          val resultsByName = response.results groupBy {
            _.name
          }
          withClue("Filter by name") { resultsByName.keys.toSet.flatten should equal(Set(Workflow1Name)) }
        }
        // Filter by multiple names
        _ <- dataAccess.queryWorkflowSummaries(
          WorkflowQueryParameters(Seq(
            WorkflowQueryKey.Name.name -> Workflow1Name, WorkflowQueryKey.Name.name -> Workflow2Name))
        ) map { case (response, _) =>
          val resultsByName = response.results groupBy {
            _.name
          }
          withClue("Filter by multiple names") { resultsByName.keys.toSet.flatten should equal(Set(Workflow1Name, Workflow2Name)) }
        }
        // Filter by workflow id
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(
          Seq(WorkflowQueryKey.Id.name -> workflow1Id.toString))) map { case (response, _) =>
          val resultsById = response.results groupBy {
            _.name
          }
          withClue("Filter by workflow ID") { resultsById.keys.toSet.flatten should equal(Set(Workflow1Name)) }
        }
        // Filter by multiple workflow ids
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(
          Seq(workflow1Id, workflow2Id).map(id => WorkflowQueryKey.Id.name -> id.toString))) map { case (response, _) =>
          val resultsById = response.results groupBy {
            _.name
          }
          withClue("Filter by multiple workflow IDs") { resultsById.keys.toSet.flatten should equal(Set(Workflow1Name, Workflow2Name)) }
        }
        // Filter by workflow id within random Ids
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(
          (randomIds :+ workflow1Id).map(id => WorkflowQueryKey.Id.name -> id.toString))) map { case (response, _) =>
          val resultsById = response.results groupBy {
            _.name
          }
          withClue("Filter by workflow ID within random IDs") { resultsById.keys.toSet.flatten should equal(Set(Workflow1Name)) }
        }
        // Filter by status
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.Status.name -> "Submitted"))) map { case (response, _) =>
          val resultsByStatus = response.results groupBy (_.status)
          withClue("Filter by status") { resultsByStatus.keys.toSet.flatten should equal(Set("Submitted")) }
        }
        // Filter by multiple statuses
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.Status.name -> "Submitted",
          WorkflowQueryKey.Status.name -> "Succeeded"))) map { case (response, _) =>
          val resultsByStatus = response.results groupBy (_.status)
          withClue("Filter by multiple statuses") { resultsByStatus.keys.toSet.flatten should equal(Set("Submitted", "Succeeded")) }
        }
        // Filter by label using AND
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.LabelAndKeyValue.name -> s"${testLabel2.key}:${testLabel2.value}"))) map { case (response, _) =>
          val resultByName = response.results groupBy (_.name)
          withClue("Filter by label using AND") { resultByName.keys.toSet.flatten should equal(Set(Workflow1Name, Workflow2Name)) }
        }
        // Filter by multiple labels using AND
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(
          Seq(testLabel2, testLabel3)
            .map(label => WorkflowQueryKey.LabelAndKeyValue.name -> s"${label.key}:${label.value}"))
        ) map { case (response, _) =>
          val resultByName = response.results groupBy (_.name)
          withClue("Filter by multiple labels using AND") { resultByName.keys.toSet.flatten should equal(Set(Workflow2Name)) }
        }
        // Filter by label using OR
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.LabelOrKeyValue.name -> s"${testLabel2.key}:${testLabel2.value}"))) map { case (response, _) =>
          val resultByName = response.results groupBy (_.name)
          withClue("Filter by label using OR") { resultByName.keys.toSet.flatten should equal(Set(Workflow1Name, Workflow2Name)) }
        }
        // Filter by multiple labels using OR
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(
          Seq(testLabel2, testLabel3)
            .map(label => WorkflowQueryKey.LabelOrKeyValue.name -> s"${label.key}:${label.value}"))
        ) map { case (response, _) =>
          val resultByName = response.results groupBy (_.name)
          withClue("Filter by multiple label using OR") { resultByName.keys.toSet.flatten should equal(Set(Workflow1Name, Workflow2Name)) }
        }
        // Filter by exclude label using AND
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.ExcludeLabelAndKeyValue.name -> s"${testLabel2.key}:${testLabel2.value}"))) map { case (response, _) =>
          val resultByName = response.results groupBy (_.name)
          withClue("Filter by exclude label using AND") {
            resultByName.keys.toSet.flatten should contain(Workflow1Name)
          }
        }
        // Filter by multiple exclude labels using AND
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(
          Seq(testLabel2, testLabel3)
            .map(label => WorkflowQueryKey.ExcludeLabelAndKeyValue.name -> s"${label.key}:${label.value}"))
        ) map { case (response, _) =>
          val resultByName = response.results groupBy (_.name)
          val ids = response.results.map(_.id)
          withClue("Filter by multiple exclude labels using AND") {
            resultByName.keys.toSet.flatten should contain(Workflow1Name)
            ids should contain(workflow1Id.toString)
            ids shouldNot contain(workflow2Id.toString)
          }
        }
        // Filter by exclude label using OR
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.ExcludeLabelOrKeyValue.name -> s"${testLabel2.key}:${testLabel2.value}"))) map { case (response, _) =>
          val resultByName = response.results groupBy (_.name)
          withClue("Filter to exclude label using OR") {
            resultByName.keys.toSet.flatten should contain(Workflow1Name)
          }
        }
        // Filter by multiple exclude labels using OR
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(
          Seq(testLabel2, testLabel3)
            .map(label => WorkflowQueryKey.ExcludeLabelOrKeyValue.name -> s"${label.key}:${label.value}"))
        ) map { case (response, _) =>
          // NOTE: On persistent databases other workflows will be returned. Just verify that our two workflows are not.
          val ids = response.results.map(_.id)
          withClue("Filter by multiple exclude labels using OR") {
            ids shouldNot contain(workflow1Id.toString)
            ids shouldNot contain(workflow2Id.toString)
          }
        }
        // Filter by start date
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.StartDate.name -> workflowQueryResult2.start.get.toUtcMilliString))) map {
          case (response, _) =>
            response.results partition {
              r => r.start.isDefined && r.start.get.compareTo(workflowQueryResult.start.get) >= 0
            } match {
              case (y, n) if y.nonEmpty && n.isEmpty => // good
              case (y, n) => fail(s"Found ${y.size} later workflows and ${n.size} earlier")
            }
        }
        // Filter by end date
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.EndDate.name -> workflowQueryResult.end.get.toUtcMilliString))) map {
          case (response, _) =>
            response.results partition {
              r => r.end.isDefined && r.end.get.compareTo(workflowQueryResult.end.get) <= 0
            } match {
              case (y, n) if y.nonEmpty && n.isEmpty => // good
              case (y, n) => fail(s"Found ${y.size} earlier workflows and ${n.size} later")
            }
        }
        // Filter by submission time
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.SubmissionTime.name -> workflowQueryResult2.submission.get.toUtcMilliString))) map {
          case (response, _) =>
            response.results partition { r =>
              r.submission.isDefined && r.submission.get.compareTo(workflowQueryResult2.submission.get) <= 0
            } match {
              case (y, n) if y.nonEmpty && n.isEmpty => // good
              case (y, n) =>
                fail(s"Found ${y.size} earlier workflows and ${n.size} later while filtering by submission timestamp")
            }
        }
        // Check for labels in query response
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.AdditionalQueryResultFields.name -> "labels"))) map {
          case (response, _) =>
            response.results filter {
              workflowQueryResult => List(workflow1Id.toString, workflow1Id.toString).contains(workflowQueryResult.id)
            } partition { r => r.labels.isDefined } match {
              case (y, n) if y.nonEmpty && n.isEmpty => //good
              case (y, n) => fail(s"Something went horribly wrong since labels were populated for ${y.size} and were missing for ${n.size} workflow(s)!")
            }
        }
      } yield ()).futureValue(Timeout(scaled(Span(30, Seconds))), Interval(scaled(Span(500, Millis))))
    }

    it should "return totalResultsCount when page and pagesize query params are specified" taggedAs DbmsTest in {
      val uniqueWorkflow3Name = s"${Workflow3Name}_${WorkflowId.randomId()}".filterNot(_ == '-')
      (for {
        _ <- baseWorkflowMetadata(uniqueWorkflow3Name)
        _ <- baseWorkflowMetadata(uniqueWorkflow3Name)
        // refresh the metadata
        _ <- dataAccess.refreshWorkflowMetadataSummaries(1000) map assertRowsProcessedAndSummarizationComplete
        //get totalResultsCount when page and pagesize are specified
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          // name being added to the query parameters so as to exclude workflows being populated by the other tests in this spec
          WorkflowQueryKey.Name.name -> uniqueWorkflow3Name,
          WorkflowQueryKey.Page.name -> "1", WorkflowQueryKey.PageSize.name -> "1"))) map { case (resp, _) =>
          resp.totalResultsCount match {
            case 2 =>
            case 1 => fail("totalResultsCount is suspiciously equal to the pageSize and not the expected total results count. Please fix!")
            case other => fail(s"totalResultsCount is expected to be 2 but is actually $other. Something has gone horribly wrong!")
          }
        }
      } yield()).futureValue
    }

    it should "revert to an prior label value" taggedAs DbmsTest in {
      def upsertLabelAndValidate(workflowId: WorkflowId, customLabelValue: String): Future[Unit] = {
        val customLabelKey = "key-1"
        val metadataEvents =
          MetadataEvent.labelsToMetadataEvents(Labels(customLabelKey -> customLabelValue), workflowId)
        for {
          _ <- dataAccess.addMetadataEvents(metadataEvents)
          _ <- dataAccess.refreshWorkflowMetadataSummaries(1000) map assertRowsProcessedAndSummarizationComplete
          _ <- dataAccess.getWorkflowLabels(workflowId).map(_.toList should contain(customLabelKey -> customLabelValue))
        } yield ()
      }

      (for {
        workflowId <- baseWorkflowMetadata("revert-label-value")
        _ <- upsertLabelAndValidate(workflowId, "value-1")
        _ <- upsertLabelAndValidate(workflowId, "")
        _ <- upsertLabelAndValidate(workflowId, "value-2")
        _ <- upsertLabelAndValidate(workflowId, "")

      } yield ()).futureValue(Timeout(scaled(Span(30, Seconds))), Interval(scaled(Span(500, Millis))))
    }

    it should "include/exclude subworkflows and check for parentWorkflowId in query response" taggedAs DbmsTest in {
      def changeParentToRunningState(id: WorkflowId): Future[Unit] = {
        val workflowKey = MetadataKey(id, jobKey = None, key = null)
        val keyAndValue = Array((WorkflowMetadataKeys.Status, WorkflowRunning.toString))

        publishMetadataEvents(workflowKey, keyAndValue)
      }

      (for {
        parentWorkflow1 <- baseWorkflowMetadata(ParentWorkflowName)
        parentWorkflow2 <- baseWorkflowMetadata(ParentWorkflow2Name)

        _ <- changeParentToRunningState(parentWorkflow1)
        // associate subworkflow 1 to parent
        _ <- subworkflowMetadata(parentWorkflow1, SubworkflowName)
        // Add in parent workflow 2
        _ <- changeParentToRunningState(parentWorkflow2)
        // associate subworkflow 2 to parent
        _ <- subworkflowMetadata(parentWorkflow2, Subworkflow2Name)
        // refresh metadata
        _ <- dataAccess.refreshWorkflowMetadataSummaries(1000)
        // include subworkflows
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(WorkflowQueryKey.IncludeSubworkflows.name -> true.toString))) map { case (resp, _) =>
          val resultByName = resp.results groupBy (_.name)
          withClue("include subworkflows:") {
            List(ParentWorkflowName, SubworkflowName, ParentWorkflow2Name, Subworkflow2Name).foreach { n => resultByName.keys.toSet.flatten should contain(n) }
          }
        }
        // exclude subworkflows - make sure we get 2 results in the right order:
        excludeSubworkflowQueryParams1 = Seq(
          WorkflowQueryKey.IncludeSubworkflows.name -> false.toString
        )
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(excludeSubworkflowQueryParams1)) map {
          case (resp, _) =>
            val resultByName = resp.results groupBy (_.name)
            withClue("exclude subworkflows (no page size):") {
              List(ParentWorkflowName, ParentWorkflow2Name).foreach { n =>
                resultByName.keys.toSet.flatten should contain(n)
              }
              List(SubworkflowName, Subworkflow2Name).foreach { n =>
                resultByName.keys.toSet.flatten should not contain n
              }
            }
        }

        // exclude subworkflows (page size 2 to make sure pagination works as expected
        excludeSubworkflowQueryParams2 = Seq(
          WorkflowQueryKey.IncludeSubworkflows.name -> false.toString,
          WorkflowQueryKey.PageSize.name -> 2.toString,
          WorkflowQueryKey.Page.name -> 1.toString
        )
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(excludeSubworkflowQueryParams2)) map {
          case (resp, _) =>
            val resultByName = resp.results groupBy (_.name)
            // This time we can do a strict equality check rather than 'contains', because we know these two should be the most recent (non-sub-) workflows:
            withClue("exclude subworkflows and assert page size:") { resultByName.keys.flatten should be(Set(ParentWorkflowName, ParentWorkflow2Name)) }
        }


        // check for parentWorkflow1 in query response
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.AdditionalQueryResultFields.name -> "parentWorkflowId",
          WorkflowQueryKey.Name.name -> SubworkflowName))) map {
          case (response, _) =>
            response.results partition { r => r.parentWorkflowId.isDefined} match {
              case (y, n) if y.nonEmpty && n.isEmpty => //good
              case (y, n) => fail(s"parentWorkflow1 should be populated for a subworkflow. It was populated correctly for ${y.size} " +
                s"and was missing in ${n.size} subworkflow(s). Something went horribly wrong!")
            }
        }
      } yield()).futureValue
    }

    it should "close the database" taggedAs DbmsTest in {
      dataAccess.metadataDatabaseInterface.close()
    }
  }
}
