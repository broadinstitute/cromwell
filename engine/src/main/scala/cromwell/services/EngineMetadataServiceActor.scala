package cromwell.services

import java.time.OffsetDateTime

import akka.actor.{Actor, ActorLogging}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.{WorkflowId, WorkflowSubmitted}
import cromwell.engine.db.DataAccess
import cromwell.engine.db.DataAccess._
import cromwell.services.EngineMetadataServiceActor._
import cromwell.services.MetadataServiceActor._
import cromwell.services.MetadataSummaryRefreshActor.{MetadataSummaryFailure, MetadataSummarySuccess, SummarizeMetadata}
import cromwell.webservice.WorkflowQueryParameters
import lenthall.config.ScalaConfig._
import spray.http.Uri

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.concurrent.duration.{Duration, FiniteDuration}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}


object EngineMetadataServiceActor {

  val MetadataSummaryRefreshInterval =
    Duration(ConfigFactory.load().getStringOr("services.MetadataService.metadata-summary-refresh-interval", "2 seconds")).asInstanceOf[FiniteDuration]

  val MetadataSummaryTimestampMinimum =
    ConfigFactory.load().getStringOption("services.MetadataService.metadata-summary-timestamp-minimum") map OffsetDateTime.parse

  // A workflow will stay in the existence cache for this many runs of the workflow summary actor before being expired out.
  val CacheExpiryCount = 5
}

// TODO: PBE: Will not be MetadataServiceActor until circular dependencies fixed.
case class EngineMetadataServiceActor(serviceConfig: Config, globalConfig: Config) extends Actor with ActorLogging {

  val dataAccess = DataAccess.globalDataAccess
  val summaryActor = context.system.actorOf(MetadataSummaryRefreshActor.props(MetadataSummaryTimestampMinimum), "metadata-summary-actor")
  self ! RefreshSummary

  // Status lookups are eventually consistent, so it's possible a db status lookup may fail for a recently submitted
  // workflow ID.  This cache records workflow IDs for which metadata has recently flowed through this actor.  If a db
  // status lookup fails, this actor consults the cache to see if the queried ID is known.  If the queried ID is known a
  // status query will return `Submitted`, otherwise the status lookup will fail.  This cache is only consulted if the
  // db status lookup fails.
  private var workflowExistenceCache: Map[WorkflowId, Int] = Map.empty

  private def queryAndRespond(query: MetadataQuery) = {
    val sndr = sender()
    dataAccess.queryMetadataEvents(query) onComplete {
      case Success(m) => sndr ! MetadataLookupResponse(query, m)
      case Failure(t) => sndr ! MetadataServiceKeyLookupFailed(query, t)
    }
  }

  private def queryStatusAndRespond(id: WorkflowId): Unit = {
    val sndr = sender()
    dataAccess.getWorkflowStatus(id) onComplete {
      case Success(Some(s)) => sndr ! StatusLookupResponse(id, s)
      case Success(None) => self ! HandleNotFound(id, sndr)
      case Failure(t) => sndr ! StatusLookupFailed(id, t)
    }
  }

  private def scheduleSummary = context.system.scheduler.scheduleOnce(MetadataSummaryRefreshInterval, self, RefreshSummary)

  private def queryWorkflowsAndRespond(uri: Uri, rawParameters: Seq[(String, String)]): Unit = {
    def queryWorkflows: Future[(WorkflowQueryResponse, Option[QueryMetadata])] = {
      for {
      // Future/Try to wrap the exception that might be thrown from WorkflowQueryParameters.apply.
        parameters <- Future.fromTry(Try(WorkflowQueryParameters(rawParameters)))
        response <- globalDataAccess.queryWorkflowSummaries(parameters)
      } yield response
    }

    val sndr = sender()

    queryWorkflows onComplete {
      case Success((response, metadata)) => sndr ! WorkflowQuerySuccess(uri, response, metadata)
      case Failure(t) => sndr ! WorkflowQueryFailure(t)
    }
  }

  def receive = {
    case action@PutMetadataAction(events) =>
      workflowExistenceCache = workflowExistenceCache ++ (events map { _.key.workflowId -> CacheExpiryCount })
      val sndr = sender()
      dataAccess.addMetadataEvents(events) onComplete {
        case Success(_) => sndr ! MetadataPutAcknowledgement(action)
        case Failure(t) =>
          val msg = MetadataPutFailed(action, t)
          log.error(t, "Sending {} failure message {}", sndr, msg)
          sndr ! msg
      }
    case GetSingleWorkflowMetadataAction(workflowId) => queryAndRespond(MetadataQuery(workflowId, None, None))
    case GetMetadataQueryAction(query@MetadataQuery(_, _, _)) => queryAndRespond(query)
    case GetStatus(workflowId) => queryStatusAndRespond(workflowId)
    case WorkflowQuery(uri, parameters) => queryWorkflowsAndRespond(uri, parameters)
    case RefreshSummary => summaryActor ! SummarizeMetadata
    case MetadataSummarySuccess =>
      // Remove expired cache entries, decrement cache counts for remaining entries.
      workflowExistenceCache = workflowExistenceCache collect { case (k, v) if v > 1 => k -> (v - 1) }
      scheduleSummary
    case MetadataSummaryFailure(t) =>
      log.error(t, "Error summarizing metadata")
      scheduleSummary
    case HandleNotFound(id, sndr) =>
      val message = if (workflowExistenceCache.contains(id)) StatusLookupResponse(id, WorkflowSubmitted) else StatusLookupNotFound(id)
      sndr ! message
  }
}
