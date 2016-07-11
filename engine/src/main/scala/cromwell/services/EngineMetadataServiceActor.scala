package cromwell.services

import java.time.OffsetDateTime
import java.util.UUID

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
      // There's a workflow existence check at the API layer.  If the request has made it this far in the system
      // then the workflow exists but it must not have generated a status yet.
      case Success(None) => sndr ! StatusLookupResponse(id, WorkflowSubmitted)
      case Failure(t) => sndr ! StatusLookupFailed(id, t)
    }
  }

  private def validateWorkflowId(validation: ValidateWorkflowIdAndExecute): Unit = {
    val possibleWorkflowId = validation.possibleWorkflowId
    val requestContext = validation.requestContext
    val callback = validation.validationCallback

    Try(UUID.fromString(possibleWorkflowId)) match {
      case Failure(t) => callback.onMalformed(possibleWorkflowId)(requestContext)
      case Success(uuid) =>
        dataAccess.workflowExistsWithId(possibleWorkflowId) onComplete {
          case Success(true) =>
            callback.onRecognized(WorkflowId(uuid))(requestContext)
          case Success(false) =>
            callback.onUnrecognized(possibleWorkflowId)(requestContext)
          case Failure(t) =>
            callback.onFailure(possibleWorkflowId, t)(requestContext)
        }
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

  private def queryWorkflowOutputsAndRespond(id: WorkflowId): Unit = {
    val replyTo = sender()
    dataAccess.queryWorkflowOutputs(id) onComplete {
      case Success(o) => replyTo ! WorkflowOutputsResponse(id, o)
      case Failure(t) => replyTo ! WorkflowOutputsFailure(id, t)
    }
  }

  private def queryLogsAndRespond(id: WorkflowId): Unit = {
    val replyTo = sender()
    dataAccess.queryLogs(id) onComplete {
      case Success(s) => replyTo ! LogsResponse(id, s)
      case Failure(t) => replyTo ! LogsFailure(id, t)
    }
  }

  def receive = {
    case action@PutMetadataAction(events) =>
      val sndr = sender()
      dataAccess.addMetadataEvents(events) onComplete {
        case Success(_) => sndr ! MetadataPutAcknowledgement(action)
        case Failure(t) =>
          val msg = MetadataPutFailed(action, t)
          log.error(t, "Sending {} failure message {}", sndr, msg)
          sndr ! msg
      }
    case v: ValidateWorkflowIdAndExecute => validateWorkflowId(v)
    case GetSingleWorkflowMetadataAction(workflowId, includeKeysOption, excludeKeysOption) =>
      queryAndRespond(MetadataQuery(workflowId, None, None, includeKeysOption, excludeKeysOption))
    case GetMetadataQueryAction(query@MetadataQuery(_, _, _, _, _)) => queryAndRespond(query)
    case GetStatus(workflowId) => queryStatusAndRespond(workflowId)
    case GetLogs(workflowId) => queryLogsAndRespond(workflowId)
    case WorkflowQuery(uri, parameters) => queryWorkflowsAndRespond(uri, parameters)
    case WorkflowOutputs(id) => queryWorkflowOutputsAndRespond(id)
    case RefreshSummary => summaryActor ! SummarizeMetadata
    case MetadataSummarySuccess => scheduleSummary
    case MetadataSummaryFailure(t) =>
      log.error(t, "Error summarizing metadata")
      scheduleSummary
  }
}
