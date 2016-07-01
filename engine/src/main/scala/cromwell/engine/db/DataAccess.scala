package cromwell.engine.db

import java.time.OffsetDateTime

import akka.actor.ActorSystem
import cromwell.core._
import cromwell.core.retry._
import cromwell.database.SqlConverters._
import cromwell.database.SqlDatabase
import cromwell.database.obj._
import cromwell.database.slick.SlickDatabase
import cromwell.engine._
import cromwell.engine.db.DataAccess.{RetryBackoff, WorkflowExecutionAndAux}
import cromwell.services.CallMetadataKeys.{ExecutionStatus => _, _}
import cromwell.services.MetadataServiceActor.{QueryMetadata, WorkflowQueryResponse}
import cromwell.services._
import cromwell.webservice.WorkflowQueryParameters
import cromwell.services.MetadataServiceActor.{QueryMetadata, WorkflowQueryResponse}
import org.slf4j.LoggerFactory
import wdl4s.Call
import cromwell.core.ExecutionIndex.IndexEnhancedIndex

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scalaz.Scalaz._
import scalaz.{NonEmptyList, Semigroup}

object DataAccess {
  lazy val log = LoggerFactory.getLogger(classOf[DataAccess])

  val globalDataAccess: DataAccess = new SlickDatabase() with DataAccess

  case class WorkflowExecutionAndAux(execution: WorkflowExecution, aux: WorkflowExecutionAux)

  val FailureEventMaxMessageLength = 1024
  val RetryBackoff = SimpleExponentialBackoff(50 millis, 1 seconds, 1D)

  private lazy val WorkflowMetadataSummarySemigroup = new Semigroup[WorkflowMetadataSummary] {
    override def append(summary1: WorkflowMetadataSummary,
                        summary2: => WorkflowMetadataSummary): WorkflowMetadataSummary = {
      // Resolve the status if both `this` and `that` have defined statuses.  This will evaluate to `None`
      // if one or both of the statuses is not defined.
      val resolvedStatus = for {
        thisStatus <- summary1.status map WorkflowState.fromString
        thatStatus <- summary2.status map WorkflowState.fromString
      } yield (thisStatus |+| thatStatus).toString

      WorkflowMetadataSummary(
        workflowUuid = summary1.workflowUuid,
        // If not both statuses are defined, take whichever is defined.
        status = resolvedStatus orElse summary1.status orElse summary2.status,
        name = summary1.name orElse summary2.name,
        startDate = summary1.startDate orElse summary2.startDate,
        endDate = summary1.endDate orElse summary2.endDate
      )
    }
  }

  private implicit class MetadatumEnhancer(val metadatum: Metadatum) extends AnyVal {
    def toSummary: WorkflowMetadataSummary = {
      val base = WorkflowMetadataSummary(metadatum.workflowUuid)
      metadatum.key match {
        case WorkflowMetadataKeys.Name => base.copy(name = metadatum.value)
        case WorkflowMetadataKeys.Status => base.copy(status = metadatum.value)
        case WorkflowMetadataKeys.StartTime =>
          base.copy(startDate = metadatum.value map OffsetDateTime.parse map { _.toSystemTimestamp })
        case WorkflowMetadataKeys.EndTime =>
          base.copy(endDate = metadatum.value map OffsetDateTime.parse map { _.toSystemTimestamp })
      }
    }
  }

  private def buildUpdatedSummary(existingSummary: Option[WorkflowMetadataSummary],
                                  metadataForUuid: Seq[Metadatum]): WorkflowMetadataSummary = {
    implicit val wmss = WorkflowMetadataSummarySemigroup

    val baseSummary = existingSummary.getOrElse(WorkflowMetadataSummary(metadataForUuid.head.workflowUuid))
    metadataForUuid.foldLeft(baseSummary) {
      case (metadataSummary, metadatum) => metadataSummary |+| metadatum.toSummary
    }
  }
}

trait DataAccess extends AutoCloseable {
  self: SqlDatabase =>

  private def withRetry[A](f: () => Future[A])(implicit actorSystem: ActorSystem): Future[A] = {
    Retry.withRetry(f, maxRetries = Option(10), backoff = RetryBackoff, isTransient = isTransient)
  }

  def addMetadataEvents(metadataEvents: Iterable[MetadataEvent])(implicit ec: ExecutionContext): Future[Unit] = {
    val metadata = metadataEvents map { metadataEvent =>
      val key = metadataEvent.key
      val workflowUuid = key.workflowId.id.toString
      val timestamp = metadataEvent.offsetDateTime.toSystemTimestamp
      val value = metadataEvent.value map { _.value }
      val valueType = metadataEvent.value map { _.valueType.typeName }
      val jobKey = key.jobKey map { jk => (jk.callFqn, jk.index, jk.attempt) }
        Metadatum(workflowUuid, key.key, jobKey.map(_._1), jobKey.flatMap(_._2), jobKey.map(_._3), value, valueType, timestamp)
    }

    addMetadata(metadata)
  }

  private def metadataToMetadataEvents(workflowId: WorkflowId)(metadata: Seq[Metadatum]): Seq[MetadataEvent] = {
    metadata map { m =>
      // If callFqn is non-null then attempt will also be non-null and there is a MetadataJobKey.
      val metadataJobKey: Option[MetadataJobKey] = for {
        callFqn <- m.callFqn
        attempt <- m.attempt
      } yield new MetadataJobKey(callFqn, m.index, attempt)

      val key = MetadataKey(workflowId, metadataJobKey, m.key)
      val value =  for {
        mValue <- m.value
        mType <- m.valueType
      } yield MetadataValue(mValue, MetadataType.fromString(mType))

      MetadataEvent(key, value, m.timestamp.toSystemOffsetDateTime)
    }
  }

  def queryMetadataEvents(query: MetadataQuery)(implicit ec: ExecutionContext): Future[Seq[MetadataEvent]] = {
    val uuid = query.workflowId.id.toString

    val futureMetadata: Future[Seq[Metadatum]] = query match {
      case MetadataQuery(_, None, None, None, None) => queryMetadataEvents(uuid)
      case MetadataQuery(_, None, Some(key), None, None) => queryMetadataEvents(uuid, key)
      case MetadataQuery(_, Some(jobKey), None, None, None) =>
        queryMetadataEvents(uuid, jobKey.callFqn, jobKey.index, jobKey.attempt)
      case MetadataQuery(_, Some(jobKey), Some(key), None, None) =>
        queryMetadataEvents(uuid, key, jobKey.callFqn, jobKey.index, jobKey.attempt)
      case MetadataQuery(_, None, None, Some(includeKeys), None) =>
        queryMetadataEventsWithWildcardKeys(uuid, includeKeys.map(_ + "%"), requireEmptyJobKey = false)
      case MetadataQuery(_, None, None, None, Some(excludeKeys)) =>
        queryMetadataEventsWithoutWildcardKeys(uuid, excludeKeys.map(_ + "%"), requireEmptyJobKey = false)
      case MetadataQuery(_, None, None, Some(includeKeys), Some(excludeKeys)) => Future.failed(
        new IllegalArgumentException(
          s"Include/Exclude keys may not be mixed: include = $includeKeys, exclude = $excludeKeys"))
      case invalidQuery => Future.failed(new IllegalArgumentException(
        s"Include/Exclude keys are only supported when querying the workflow, not when querying calls: $invalidQuery"))
    }

    futureMetadata map metadataToMetadataEvents(query.workflowId)
  }

  def queryWorkflowOutputs(id: WorkflowId)
                          (implicit ec: ExecutionContext): Future[Seq[MetadataEvent]] = {
    val uuid = id.id.toString
    queryMetadataEventsWithWildcardKeys(uuid, s"${WorkflowMetadataKeys.Outputs}:%".wrapNel, requireEmptyJobKey = true).
      map(metadataToMetadataEvents(id))
  }

  def queryLogs(id: WorkflowId)
               (implicit ec: ExecutionContext): Future[Seq[MetadataEvent]] = {

    val uuid = id.id.toString
    val keys = NonEmptyList(Stdout, Stderr, BackendLogsPrefix + ":%")
    queryMetadataEventsWithWildcardKeys(id.id.toString, keys, requireEmptyJobKey = false) map metadataToMetadataEvents(id)
  }

  def refreshWorkflowMetadataSummaries(startMetadataId: Long, startMetadataDateTime: Option[OffsetDateTime])
                                      (implicit ec: ExecutionContext): Future[Long] = {
    self.refreshMetadataSummaries(startMetadataId, startMetadataDateTime map { _.toSystemTimestamp },
      DataAccess.buildUpdatedSummary)
  }

  def getWorkflowStatus(id: WorkflowId)
                       (implicit ec: ExecutionContext): Future[Option[WorkflowState]] = {
    self.getStatus(id.toString) map { _ map WorkflowState.fromString }
  }

  def queryWorkflowSummaries(queryParameters: WorkflowQueryParameters)
                            (implicit ec: ExecutionContext): Future[(WorkflowQueryResponse, Option[QueryMetadata])] = {
    val workflowSummaries = queryWorkflowSummaries(
      queryParameters.statuses, queryParameters.names, queryParameters.ids.map(_.toString),
      queryParameters.startDate.map(_.toSystemTimestamp), queryParameters.endDate.map(_.toSystemTimestamp),
      queryParameters.page, queryParameters.pageSize)

    val workflowSummaryCount = countWorkflowSummaries(
      queryParameters.statuses, queryParameters.names, queryParameters.ids.map(_.toString),
      queryParameters.startDate.map(_.toSystemTimestamp), queryParameters.endDate.map(_.toSystemTimestamp))

    workflowSummaryCount flatMap { count =>
      workflowSummaries map { workflows =>
        (WorkflowQueryResponse(workflows.toSeq map { workflow =>
          MetadataServiceActor.WorkflowQueryResult(
            id = workflow.workflowUuid,
            name = workflow.name,
            status = workflow.status,
            start = workflow.startDate map { _.toSystemOffsetDateTime },
            end = workflow.endDate map { _.toSystemOffsetDateTime })
        }),
          //only return metadata if page is defined
          queryParameters.page map { _ => QueryMetadata(queryParameters.page, queryParameters.pageSize, Option(count)) })
      }
    }
  }

  def getExecutionInfoByKey(workflowId: WorkflowId, call: Call, attempt: Int, key: String)
                           (implicit ec: ExecutionContext): Future[Option[Option[String]]] = {
    getExecutionInfoByKey(workflowId.toString, call.fullyQualifiedName, attempt, key)
  }

  /**
    * TODO: the interface for retrying SQL commands that might fail because
    * of transient reasons (e.g. upsertExecutionInfo and upsertRuntimeAttributes)
    * could be made better.  Also it'd be nice if it were transparently
    * turned on/off for all methods in dataAccess.
    *
    * https://github.com/broadinstitute/cromwell/issues/693
    */
  def upsertExecutionInfo(workflowId: WorkflowId,
                          callKey: JobKey,
                          keyValues: Map[String, Option[String]],
                          actorSystem: ActorSystem): Future[Unit] = {
    implicit val system = actorSystem
    implicit val ec = actorSystem.dispatcher
    withRetry(() => upsertExecutionInfo(workflowId.toString, callKey.scope.fullyQualifiedName, callKey.index.fromIndex,
      callKey.attempt, keyValues))
  }
}
