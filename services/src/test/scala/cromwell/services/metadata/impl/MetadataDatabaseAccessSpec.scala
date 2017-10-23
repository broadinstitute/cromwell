package cromwell.services.metadata.impl

import java.time.OffsetDateTime

import com.typesafe.config.ConfigFactory
import cromwell.core.Tags.DbmsTest
import cromwell.core._
import cromwell.core.labels.Label
import cromwell.database.slick.MetadataSlickDatabase
import cromwell.services.MetadataServicesStore
import cromwell.services.ServicesStore.EnhancedSqlDatabase
import cromwell.services.metadata._
import org.scalatest.concurrent.PatienceConfiguration.{Interval, Timeout}
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import org.specs2.mock.Mockito

import scala.concurrent.{ExecutionContext, Future}
object MetadataDatabaseAccessSpec {
  val AllowFalse = Seq(QueryParameter("allow", "false"))
  val AllowTrue = Seq(QueryParameter("allow", "true"))

  val Workflow1Name = "test1"
  val Workflow2Name = "test2"
}

class MetadataDatabaseAccessSpec extends FlatSpec with Matchers with ScalaFutures with BeforeAndAfterAll with Mockito {
  import MetadataDatabaseAccessSpec._

  "MetadataDatabaseAccess (hsqldb)" should behave like testWith("database")

  "MetadataDatabaseAccess (mysql)" should behave like testWith("database-test-mysql")

  implicit val ec = ExecutionContext.global

  implicit val defaultPatience = PatienceConfig(scaled(Span(5, Seconds)), scaled(Span(100, Millis)))

  def testWith(configPath: String): Unit = {
    lazy val dataAccess = new MetadataDatabaseAccess with MetadataServicesStore {
      override val metadataDatabaseInterface =
        new MetadataSlickDatabase(ConfigFactory.load.getConfig(configPath))
          .initialized(MetadataServicesStore.MetadataLiquibaseSettings)
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
        (WorkflowMetadataKeys.StartTime, OffsetDateTime.now.toString),
        (WorkflowMetadataKeys.Status, WorkflowSubmitted.toString),
        (WorkflowMetadataKeys.Name, name)
      ) ++ labelMetadata

      publishMetadataEvents(workflowKey, keyAndValue(name)).map(_ => workflowId)
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
        response <- dataAccess.queryMetadataEvents(MetadataQuery(workflow1Id, None, Option(WorkflowMetadataKeys.WorkflowRoot), None, None, expandSubWorkflows = false))
        _ = response shouldBe expected
      } yield()).futureValue
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
          (WorkflowMetadataKeys.EndTime, OffsetDateTime.now.toString))

        publishMetadataEvents(workflowKey, keyAndValue)
      }

      (for {
        workflow1Id <- baseWorkflowMetadata(Workflow1Name, Set(testLabel1, testLabel2))
        _ <- succeededWorkflowMetadata(workflow1Id)
        // Put a bit of space between the two workflows
        _ = Thread.sleep(50)
        workflow2Id <- baseWorkflowMetadata(Workflow2Name, Set(testLabel2, testLabel3))

        // refresh the metadata
        _ <- dataAccess.refreshWorkflowMetadataSummaries() map { max =>
          max should be > 0L
        }

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
          resultsByName.keys.toSet.flatten should equal(Set(Workflow1Name))
        }
        // Filter by multiple names
        _ <- dataAccess.queryWorkflowSummaries(
          WorkflowQueryParameters(Seq(
            WorkflowQueryKey.Name.name -> Workflow1Name, WorkflowQueryKey.Name.name -> Workflow2Name))
        ) map { case (response, _) =>
          val resultsByName = response.results groupBy {
            _.name
          }
          resultsByName.keys.toSet.flatten should equal(Set(Workflow1Name, Workflow2Name))
        }
        // Filter by workflow id
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(
          Seq(WorkflowQueryKey.Id.name -> workflow1Id.toString))) map { case (response, _) =>
          val resultsById = response.results groupBy {
            _.name
          }
          resultsById.keys.toSet.flatten should equal(Set(Workflow1Name))
        }
        // Filter by multiple workflow ids
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(
          Seq(workflow1Id, workflow2Id).map(id => WorkflowQueryKey.Id.name -> id.toString))) map { case (response, _) =>
          val resultsById = response.results groupBy {
            _.name
          }
          resultsById.keys.toSet.flatten should equal(Set(Workflow1Name, Workflow2Name))
        }
        // Filter by workflow id within random Ids
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(
          (randomIds :+ workflow1Id).map(id => WorkflowQueryKey.Id.name -> id.toString))) map { case (response, _) =>
          val resultsById = response.results groupBy {
            _.name
          }
          resultsById.keys.toSet.flatten should equal(Set(Workflow1Name))
        }
        // Filter by status
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.Status.name -> "Submitted"))) map { case (response, _) =>
          val resultsByStatus = response.results groupBy (_.status)
          resultsByStatus.keys.toSet.flatten should equal(Set("Submitted"))
        }
        // Filter by multiple statuses
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.Status.name -> "Submitted",
          WorkflowQueryKey.Status.name -> "Succeeded"))) map { case (response, _) =>
          val resultsByStatus = response.results groupBy (_.status)
          resultsByStatus.keys.toSet.flatten should equal(Set("Submitted", "Succeeded"))
        }
        // Filter by label
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.LabelKeyValue.name -> s"${testLabel2.key}:${testLabel2.value}"))) map { case (response, _) =>
          val resultByName = response.results groupBy (_.name)
          resultByName.keys.toSet.flatten should equal(Set(Workflow1Name, Workflow2Name))
        }
        // Filter by multiple labels
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(
          Seq(testLabel2, testLabel3)
            .map(label => WorkflowQueryKey.LabelKeyValue.name -> s"${label.key}:${label.value}"))
        ) map { case (response, _) =>
          val resultByName = response.results groupBy (_.name)
          resultByName.keys.toSet.flatten should equal(Set(Workflow2Name))
        }
        // Filter by start date
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.StartDate.name -> workflowQueryResult2.start.get.toString))) map { case (response, _) =>
          response.results partition { r => r.start.isDefined && r.start.get.compareTo(workflowQueryResult.start.get) >= 0 } match {
            case (y, n) if y.nonEmpty && n.isEmpty => // good
            case (y, n) => fail(s"Found ${y.size} later workflows and ${n.size} earlier")
          }
        }
        // Filter by end date
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.EndDate.name -> workflowQueryResult.end.get.toString))) map { case (response, _) =>
          response.results partition { r => r.end.isDefined && r.end.get.compareTo(workflowQueryResult.end.get) <= 0 } match {
            case (y, n) if y.nonEmpty && n.isEmpty => // good
            case (y, n) => fail(s"Found ${y.size} earlier workflows and ${n.size} later")
          }
        }
      } yield ()).futureValue(Timeout(scaled(Span(30, Seconds))), Interval(scaled(Span(500, Millis))))
    }

    it should "close the database" taggedAs DbmsTest in {
      dataAccess.metadataDatabaseInterface.close()
    }
  }
}
