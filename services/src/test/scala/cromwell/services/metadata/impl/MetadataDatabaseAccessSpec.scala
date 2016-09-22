package cromwell.services.metadata.impl

import java.time.OffsetDateTime

import com.typesafe.config.ConfigFactory
import cromwell.core.Tags.DbmsTest
import cromwell.core._
import cromwell.database.slick.SlickDatabase
import cromwell.services.ServicesStore
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

    import ServicesStore.EnhancedSqlDatabase

    lazy val dataAccess: MetadataDatabaseAccess = new MetadataDatabaseAccess with ServicesStore {
      override val databaseInterface = new SlickDatabase(ConfigFactory.load.getConfig(configPath)).initialized
    }

    def publishMetadataEvents(baseKey: MetadataKey, keyValues: Array[(String, String)]): Future[Unit] = {
      val events = keyValues map { case (k, v) =>
        MetadataEvent(baseKey.copy(key = k), MetadataValue(v))
      }
      dataAccess.addMetadataEvents(events)
    }

    def baseWorkflowMetadata(name: String): Future[WorkflowId] = {
      val workflowId = WorkflowId.randomId()
      val workflowKey = MetadataKey(workflowId, jobKey = None, key = null)
      def keyAndValue(name: String) = Array(
        (WorkflowMetadataKeys.StartTime, OffsetDateTime.now.toString),
        (WorkflowMetadataKeys.Status, WorkflowSubmitted.toString),
        (WorkflowMetadataKeys.Name, name)
      )
      publishMetadataEvents(workflowKey, keyAndValue(name)).map(_ => workflowId)
    }

    it should "return pagination metadata only when page and pagesize query params are specified" taggedAs DbmsTest in {
      (for {
        workflow1Id <- baseWorkflowMetadata(Workflow1Name)
        //get metadata when page and pagesize are specified
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.Page.name -> "1", WorkflowQueryKey.PageSize.name -> "50"))) map { case (response, meta) =>
          meta match {
            case Some(metadata) =>
            case None => fail("Should have metadata when page and pagesize are specified.")
          }
        }
        //don't get metadata when page and pagesize are not specified
        _ <- dataAccess.queryWorkflowSummaries(
          WorkflowQueryParameters(Seq())) map { case(response, meta) =>
          meta match {
            case Some(metadata) => fail("Should not have metadata when page and pagesize are not specified")
            case None =>
          }
        }
      } yield()).futureValue
    }

    it should "create and query a workflow" taggedAs DbmsTest in {

      val randomIds = Seq.fill(10)(WorkflowId.randomId().toString)

      def succeededWorkflowMetadata(id: WorkflowId): Future[Unit] = {
        val workflowKey = MetadataKey(id, jobKey = None, key = null)
        val keyAndValue = Array(
          (WorkflowMetadataKeys.Status, WorkflowRunning.toString),
          (WorkflowMetadataKeys.Status, WorkflowSucceeded.toString),
          (WorkflowMetadataKeys.EndTime, OffsetDateTime.now.toString))

        publishMetadataEvents(workflowKey, keyAndValue)
      }

      (for {
        workflow1Id <- baseWorkflowMetadata(Workflow1Name)
        _ <- succeededWorkflowMetadata(workflow1Id)
        // Put a bit of space between the two workflows
        _ = Thread.sleep(50)
        workflow2Id <- baseWorkflowMetadata(Workflow2Name)

        // refresh the metadata
        _ <- dataAccess.refreshWorkflowMetadataSummaries() map { max =>
          max should be > 0L
        }

        // Query with no filters
        (workflowQueryResult, workflowQueryResult2) <-
        dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq.empty)) map { case (response, meta) =>
          val result = response.results find { r => r.name.contains(Workflow1Name) && r.end.isDefined } getOrElse
            fail(s"$Workflow1Name with an end not found in ${response.results}")
          val result2 = response.results find {
            _.name.contains(Workflow2Name)
          } getOrElse fail(s"$Workflow2Name not found in ${response.results}")
          (result, result2)
        }
        // Filter by name
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(WorkflowQueryKey.Name.name -> Workflow1Name))) map { case (response, meta) =>
          val resultsByName = response.results groupBy {
            _.name
          }
          resultsByName.keys.toSet.flatten should equal(Set(Workflow1Name))
        }
        // Filter by multiple names
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(WorkflowQueryKey.Name.name -> Workflow1Name, WorkflowQueryKey.Name.name -> Workflow2Name))) map { case (response, meta) =>
          val resultsByName = response.results groupBy {
            _.name
          }
          resultsByName.keys.toSet.flatten should equal(Set(Workflow1Name, Workflow2Name))
        }
        // Filter by workflow id
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(
          Seq(WorkflowQueryKey.Id.name -> workflow1Id.toString))) map { case (response, meta) =>
          val resultsById = response.results groupBy {
            _.name
          }
          resultsById.keys.toSet.flatten should equal(Set(Workflow1Name))
        }
        // Filter by multiple workflow ids
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(
          Seq(workflow1Id, workflow2Id).map(id => WorkflowQueryKey.Id.name -> id.toString))) map { case (response, meta) =>
          val resultsById = response.results groupBy {
            _.name
          }
          resultsById.keys.toSet.flatten should equal(Set(Workflow1Name, Workflow2Name))
        }
        // Filter by workflow id within random Ids
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(
          (randomIds :+ workflow1Id).map(id => WorkflowQueryKey.Id.name -> id.toString))) map { case (response, meta) =>
          val resultsById = response.results groupBy {
            _.name
          }
          resultsById.keys.toSet.flatten should equal(Set(Workflow1Name))
        }
        // Filter by status
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(WorkflowQueryKey.Status.name -> "Submitted"))) map { case (response, meta) =>
          val resultsByStatus = response.results groupBy (_.status)
          resultsByStatus.keys.toSet.flatten should equal(Set("Submitted"))
        }
        // Filter by multiple statuses
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(WorkflowQueryKey.Status.name -> "Submitted", WorkflowQueryKey.Status.name -> "Succeeded"))) map { case (response, meta) =>
          val resultsByStatus = response.results groupBy (_.status)
          resultsByStatus.keys.toSet.flatten should equal(Set("Submitted", "Succeeded"))
        }
        // Filter by start date
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.StartDate.name -> workflowQueryResult2.start.get.toString))) map { case (response, meta) =>
          response.results partition { r => r.start.isDefined && r.start.get.compareTo(workflowQueryResult.start.get) >= 0 } match {
            case (y, n) if y.nonEmpty && n.isEmpty => // good
            case (y, n) => fail(s"Found ${y.size} later workflows and ${n.size} earlier")
          }
        }
        // Filter by end date
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.EndDate.name -> workflowQueryResult.end.get.toString))) map { case (response, meta) =>
          response.results partition { r => r.end.isDefined && r.end.get.compareTo(workflowQueryResult.end.get) <= 0 } match {
            case (y, n) if y.nonEmpty && n.isEmpty => // good
            case (y, n) => fail(s"Found ${y.size} earlier workflows and ${n.size} later")
          }
        }
      } yield ()).futureValue(Timeout(scaled(Span(30, Seconds))), Interval(scaled(Span(500, Millis))))
    }

    it should "close the database" taggedAs DbmsTest in {
      dataAccess.closeDatabaseInterface()
    }
  }
}
