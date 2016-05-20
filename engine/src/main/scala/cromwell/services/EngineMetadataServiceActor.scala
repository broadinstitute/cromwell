package cromwell.services

import java.time.OffsetDateTime

import akka.actor.{Actor, ActorLogging}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.WorkflowId
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
}

// TODO: PBE: Will not be MetadataServiceActor until circular dependencies fixed.
case class EngineMetadataServiceActor(serviceConfig: Config, globalConfig: Config) extends Actor with ActorLogging {

  val dataAccess = DataAccess.globalDataAccess
  val summaryActor = context.system.actorOf(MetadataSummaryRefreshActor.props(MetadataSummaryTimestampMinimum), "metadata-summary-actor")
  self ! RefreshSummary

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
      case Success(None) => sndr ! StatusLookupNotFound(id)
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
    case action@PutMetadataAction(event) =>
      val sndr = sender()
      dataAccess.addMetadataEvent(event) onComplete {
        case Success(_) => sndr ! MetadataPutAcknowledgement(action)
        case Failure(t) =>
          val msg = MetadataPutFailed(action, t)
          log.error(t, "Sending {} failure message {}", sndr, msg)
          sndr ! msg
      }
    case GetAllMetadataAction(workflowId) => queryAndRespond(MetadataQuery(workflowId, None, None))
    case GetMetadataQueryAction(query@MetadataQuery(_, _, _)) => queryAndRespond(query)
    case GetStatus(workflowId) => queryStatusAndRespond(workflowId)
    case WorkflowQuery(uri, parameters) => queryWorkflowsAndRespond(uri, parameters)
    case RefreshSummary => summaryActor ! SummarizeMetadata
    case MetadataSummarySuccess => scheduleSummary
    case MetadataSummaryFailure(t) =>
      log.error(t, "Error summarizing metadata")
      scheduleSummary
  }
}
