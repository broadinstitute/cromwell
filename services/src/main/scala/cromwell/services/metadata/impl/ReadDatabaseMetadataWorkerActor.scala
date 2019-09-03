package cromwell.services.metadata.impl

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
  def props(metadataReadTimeout: Duration) = Props(new ReadDatabaseMetadataWorkerActor(metadataReadTimeout)).withDispatcher(ServiceDispatcher)
}

class ReadDatabaseMetadataWorkerActor(metadataReadTimeout: Duration) extends Actor with ActorLogging with MetadataDatabaseAccess with MetadataServicesStore {

  implicit val ec = context.dispatcher

  def receive = {
    case GetMetadataAction(query@MetadataQuery(_, _, _, _, _, _)) => evaluateRespondAndStop(sender(), getMetadata(query))
    case GetStatus(workflowId) => evaluateRespondAndStop(sender(), getStatus(workflowId))
    case GetLabels(workflowId) => evaluateRespondAndStop(sender(), queryLabelsAndRespond(workflowId))
    case GetLogs(workflowId) => evaluateRespondAndStop(sender(), queryLogsAndRespond(workflowId))
    case query: QueryForWorkflowsMatchingParameters => evaluateRespondAndStop(sender(), queryWorkflowsAndRespond(query.parameters))
    case WorkflowOutputs(id) => evaluateRespondAndStop(sender(), queryWorkflowOutputsAndRespond(id))
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

  private def getMetadata(query: MetadataQuery): Future[MetadataServiceResponse] = {

    queryMetadataEvents(query, metadataReadTimeout) map {
      m => MetadataLookupResponse(query, m)
    } recover {
      case t => MetadataServiceKeyLookupFailed(query, t)
    }
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
