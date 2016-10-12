package cromwell.services.metadata.impl

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.core.Dispatcher.ApiDispatcher
import cromwell.core.{WorkflowId, WorkflowSubmitted}
import cromwell.services.SingletonServicesStore
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata.{CallMetadataKeys, MetadataQuery, WorkflowQueryParameters}

import scala.concurrent.Future
import scala.util.{Failure, Success, Try}

object ReadMetadataActor {
  def props() = Props(new ReadMetadataActor()).withDispatcher(ApiDispatcher)
}

class ReadMetadataActor extends Actor with ActorLogging with MetadataDatabaseAccess with SingletonServicesStore {

  implicit val ec = context.dispatcher

  def receive = {
    case GetSingleWorkflowMetadataAction(workflowId, includeKeysOption, excludeKeysOption, expandSubWorkflows) =>
      val includeKeys = if (expandSubWorkflows) {
        includeKeysOption map { _.::(CallMetadataKeys.SubWorkflowId) }
      } else includeKeysOption
      queryAndRespond(MetadataQuery(workflowId, None, None, includeKeys, excludeKeysOption, expandSubWorkflows))
    case GetMetadataQueryAction(query@MetadataQuery(_, _, _, _, _, _)) => queryAndRespond(query)
    case GetStatus(workflowId) => queryStatusAndRespond(workflowId)
    case GetLogs(workflowId) => queryLogsAndRespond(workflowId)
    case query: WorkflowQuery[_] => queryWorkflowsAndRespond(query.uri, query.parameters)
    case WorkflowOutputs(id) => queryWorkflowOutputsAndRespond(id)
  }

  private def queryAndRespond(query: MetadataQuery) = {
    val sndr = sender()
    queryMetadataEvents(query) onComplete {
      case Success(m) => sndr ! MetadataLookupResponse(query, m)
      case Failure(t) => sndr ! MetadataServiceKeyLookupFailed(query, t)
    }
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

  private def queryWorkflowsAndRespond[A](uri: A, rawParameters: Seq[(String, String)]): Unit = {
    def queryWorkflows: Future[(WorkflowQueryResponse, Option[QueryMetadata])] = {
      for {
      // Future/Try to wrap the exception that might be thrown from WorkflowQueryParameters.apply.
        parameters <- Future.fromTry(Try(WorkflowQueryParameters(rawParameters)))
        response <- queryWorkflowSummaries(parameters)
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
    queryWorkflowOutputs(id) onComplete {
      case Success(o) => replyTo ! WorkflowOutputsResponse(id, o)
      case Failure(t) => replyTo ! WorkflowOutputsFailure(id, t)
    }
  }

  private def queryLogsAndRespond(id: WorkflowId): Unit = {
    val replyTo = sender()
    queryLogs(id) onComplete {
      case Success(s) => replyTo ! LogsResponse(id, s)
      case Failure(t) => replyTo ! LogsFailure(id, t)
    }
  }

}
