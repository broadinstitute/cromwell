package cromwell.services.metadata.impl

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.core.Dispatcher.ApiDispatcher
import cromwell.core.{WorkflowId, WorkflowSubmitted}
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata.WorkflowQueryParameters

import scala.concurrent.Future
import scala.util.{Failure, Success, Try}

object ReadMetadataActor {
  def props() = Props(new ReadMetadataActor()).withDispatcher(ApiDispatcher)
}

class ReadMetadataActor extends Actor with ActorLogging with MetadataDatabaseAccess with MetadataServicesStore {

  implicit val ec = context.dispatcher

  def receive = {
    case GetStatus(workflowId) => queryStatusAndRespond(workflowId)
    case GetLabels(workflowId) => queryLabelsAndRespond(workflowId)
    case query: WorkflowQuery => queryWorkflowsAndRespond(query.parameters)
  }

  private def queryStatusAndRespond(id: WorkflowId): Unit = {
    val sndr = sender()
    getWorkflowStatus(id) onComplete {
      case Success(Some(s)) => sndr ! StatusLookupResponse(id, s)
      // There's a workflow existence check at the API layer.  If the request has made it this far in the system
      // then the workflow exists but it must not have generated a status yet.
      case Success(None) => sndr ! StatusLookupResponse(id, WorkflowSubmitted)
      case Failure(t) => sndr ! StatusLookupFailed(id, t)
    }
  }

  private def queryLabelsAndRespond(id: WorkflowId): Unit = {
    val sndr = sender()
    getWorkflowLabels(id) onComplete {
      case Success(ls) => sndr ! LabelLookupResponse(id, ls)
      case Failure(t) => sndr ! LabelLookupFailed(id, t)
    }
  }

  private def queryWorkflowsAndRespond(rawParameters: Seq[(String, String)]): Unit = {
    def queryWorkflows: Future[(WorkflowQueryResponse, Option[QueryMetadata])] = {
      for {
      // Future/Try to wrap the exception that might be thrown from WorkflowQueryParameters.apply.
        parameters <- Future.fromTry(Try(WorkflowQueryParameters(rawParameters)))
        response <- queryWorkflowSummaries(parameters)
      } yield response
    }

    val sndr = sender()

    queryWorkflows onComplete {
      case Success((response, metadata)) => sndr ! WorkflowQuerySuccess(response, metadata)
      case Failure(t) => sndr ! WorkflowQueryFailure(t)
    }
  }
}
