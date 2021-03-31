package cromwell.services.metadata.impl

import java.sql.SQLTimeoutException

import akka.actor.{Actor, ActorLogging, ActorRef, PoisonPill, Props}
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.core.{WorkflowId, WorkflowSubmitted}
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata.{MetadataQuery, WorkflowQueryParameters}

import scala.concurrent.Future
import scala.concurrent.duration.Duration
import scala.util.Try

object ReadDatabaseMetadataWorkerActor {
  def props(metadataReadTimeout: Duration, metadataReadRowNumberSafetyThreshold: Int) =
    Props(new ReadDatabaseMetadataWorkerActor(metadataReadTimeout, metadataReadRowNumberSafetyThreshold)).withDispatcher(ServiceDispatcher)
}

class ReadDatabaseMetadataWorkerActor(metadataReadTimeout: Duration, metadataReadRowNumberSafetyThreshold: Int)
  extends Actor
    with ActorLogging
    with MetadataDatabaseAccess
    with MetadataServicesStore {

  implicit val ec = context.dispatcher

  def receive = {
    case GetMetadataAction(query: MetadataQuery, checkTotalMetadataRowNumberBeforeQuerying: Boolean) =>
      evaluateRespondAndStop(sender(), getMetadata(query, checkTotalMetadataRowNumberBeforeQuerying))
    case GetMetadataStreamAction(workflowId, fetchSize) =>
      evaluateRespondAndStop(sender(), Future.fromTry(getMetadataStream(workflowId, fetchSize)))
    case GetStatus(workflowId) => evaluateRespondAndStop(sender(), getStatus(workflowId))
    case GetLabels(workflowId) => evaluateRespondAndStop(sender(), queryLabelsAndRespond(workflowId))
    case GetRootAndSubworkflowLabels(rootWorkflowId: WorkflowId) => evaluateRespondAndStop(sender(), queryRootAndSubworkflowLabelsAndRespond(rootWorkflowId))
    case GetLogs(workflowId) => evaluateRespondAndStop(sender(), queryLogsAndRespond(workflowId))
    case QueryForWorkflowsMatchingParameters(parameters) => evaluateRespondAndStop(sender(), queryWorkflowsAndRespond(parameters))
    case WorkflowOutputs(id) => evaluateRespondAndStop(sender(), queryWorkflowOutputsAndRespond(id))
    case unexpected => log.warning(s"Programmer Error! Unexpected message received by ${getClass.getSimpleName}: $unexpected")
  }

  private def evaluateRespondAndStop(sndr: ActorRef, f: Future[Any]) = {
    f map { result =>
      sndr ! result
    } andThen {
      case _ => self ! PoisonPill
    } recover {
      case t => log.error(t, s"Programmer Error! Unexpected error fall-through to 'evaluateRespondAndStop in ${getClass.getSimpleName}'")
    }
    ()
  }

  private def getMetadata(query: MetadataQuery, checkResultSizeBeforeQuerying: Boolean): Future[MetadataServiceResponse] = {
    if (checkResultSizeBeforeQuerying) {
      getMetadataReadRowCount(query, metadataReadTimeout) flatMap { count =>
        if (count > metadataReadRowNumberSafetyThreshold) {
          Future.successful(MetadataLookupFailedTooLargeResponse(query, count))
        } else {
          queryMetadata(query)
        }
      } recoverWith {
        case _: SQLTimeoutException => Future.successful(MetadataLookupFailedTimeoutResponse(query))
        case t => Future.successful(MetadataServiceKeyLookupFailed(query, t))
      }
    } else {
      queryMetadata(query)
    }
  }

  private def getMetadataStream(workflowId: WorkflowId, fetchSize: Int): Try[MetadataServiceResponse] = {
    metadataEventsStream(workflowId, fetchSize) map {
      s => MetadataLookupStreamSuccess(workflowId, s)
    } recover {
      case t => MetadataLookupStreamFailed(workflowId, t)
    }
  }

  private def queryMetadata(query: MetadataQuery): Future[MetadataServiceResponse] =
    queryMetadataEvents(query, metadataReadTimeout) map {
      m => MetadataLookupResponse(query, m)
    } recover {
      case _: SQLTimeoutException => MetadataLookupFailedTimeoutResponse(query)
      case t => MetadataServiceKeyLookupFailed(query, t)
    }

  private def getStatus(id: WorkflowId): Future[MetadataServiceResponse] = {

    getWorkflowStatus(id) map {
      case Some(s) => StatusLookupResponse(id, s)
      // There's a workflow existence check at the API layer.  If the request has made it this far in the system
      // then the workflow exists but it must not have generated a status yet.
      case None => StatusLookupResponse(id, WorkflowSubmitted)
    } recover {
      case t => StatusLookupFailed(id, t)
    }
  }

  private def queryLabelsAndRespond(id: WorkflowId): Future[MetadataServiceResponse] = {

    getWorkflowLabels(id) map {
      ls => LabelLookupResponse(id, ls)
    } recover {
      case t => LabelLookupFailed(id, t)
    }
  }

  private def queryRootAndSubworkflowLabelsAndRespond(rootWorkflowId: WorkflowId): Future[MetadataServiceResponse] = {

    getRootAndSubworkflowLabels(rootWorkflowId) map { labels =>
      RootAndSubworkflowLabelsLookupResponse(rootWorkflowId, labels)
    } recover {
      case t => RootAndSubworkflowLabelsLookupFailed(rootWorkflowId, t)
    }
  }

  private def queryWorkflowsAndRespond(rawParameters: Seq[(String, String)]): Future[MetadataServiceResponse] = {
    def queryWorkflows: Future[(WorkflowQueryResponse, Option[QueryMetadata])] = {
      for {
      // Future/Try to wrap the exception that might be thrown from WorkflowQueryParameters.apply.
        parameters <- Future.fromTry(Try(WorkflowQueryParameters(rawParameters)))
        response <- queryWorkflowSummaries(parameters)
      } yield response
    }

    queryWorkflows map {
      case (response, metadata) => WorkflowQuerySuccess(response, metadata)
    } recover {
      case t => WorkflowQueryFailure(t)
    }
  }

  private def queryWorkflowOutputsAndRespond(id: WorkflowId): Future[MetadataServiceResponse] = {
    queryWorkflowOutputs(id, metadataReadTimeout) map {
      o => WorkflowOutputsResponse(id, o)
    } recover {
      case t => WorkflowOutputsFailure(id, t)
    }
  }

  private def queryLogsAndRespond(id: WorkflowId): Future[MetadataServiceResponse] = {
    queryLogs(id, metadataReadTimeout) map {
      s => LogsResponse(id, s)
    } recover {
      case t => LogsFailure(id, t)
    }
  }

}
