package cromwell.engine.workflow.lifecycle.finalization

import akka.actor.{Actor, ActorLogging, ActorRef, ActorSystem, Props}
import akka.event.LoggingReceive
import akka.http.scaladsl.Http
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.marshalling.Marshal
import akka.http.scaladsl.model._
import akka.util.ByteString
import cromwell.backend.BackendLifecycleActor.BackendWorkflowLifecycleActorResponse
import cromwell.backend.BackendWorkflowFinalizationActor.{FinalizationFailed, FinalizationSuccess, Finalize}
import cromwell.core.Dispatcher.IoDispatcher
import cromwell.core.retry.Retry.withRetry
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.core.{CallOutputs, WorkflowId, WorkflowState}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.finalization.WorkflowCallbackJsonSupport._

import scala.concurrent.ExecutionContextExecutor
//import cromwell.services.ServiceRegistryActor

import java.net.URI
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

/**
  * The WorkflowCallbackActor is responsible for sending a message on workflow completion to a configured endpoint.
  * This allows for users to build automation around workflow completion without needing to poll Cromwell for
  * workflow status.
  */
object WorkflowCallbackActor {
  def props(workflowId: WorkflowId,
           // serviceRegistryActor: ServiceRegistryActor,
            workflowDescriptor: EngineWorkflowDescriptor,
            workflowOutputs: CallOutputs,
            lastWorkflowState: WorkflowState,
            callbackUri: URI,
            retryBackoff: SimpleExponentialBackoff
           ) = Props(
    new WorkflowCallbackActor(
      workflowId,
      /*serviceRegistryActor,*/
      workflowDescriptor,
      workflowOutputs,
      lastWorkflowState,
      callbackUri,
      retryBackoff)
  ).withDispatcher(IoDispatcher)
}

class WorkflowCallbackActor(workflowId: WorkflowId,
                            //serviceRegistryActor: ServiceRegistryActor,
                            val workflowDescriptor: EngineWorkflowDescriptor,
                            workflowOutputs: CallOutputs,
                            lastWorkflowState: WorkflowState,
                            callbackUri: URI,
                            retryBackoff: SimpleExponentialBackoff)
  extends Actor with ActorLogging {

  implicit val ec: ExecutionContextExecutor = context.dispatcher
  implicit val system: ActorSystem = context.system

  lazy private val callbackMessage = CallbackMessage(
    workflowId.toString, lastWorkflowState.toString, workflowOutputs.outputs.map(entry => (entry._1.name, entry._2))
  )

  override def receive = LoggingReceive {
    case Finalize => performActionThenRespond(
      performCallback().map( _ => FinalizationSuccess),
      FinalizationFailed
    )
  }

  private def performActionThenRespond(operation: => Future[BackendWorkflowLifecycleActorResponse],
                                       onFailure: (Throwable) => BackendWorkflowLifecycleActorResponse)
                                      (implicit ec: ExecutionContext) = {
    val respondTo: ActorRef = sender()
    operation onComplete {
      case Success(r) => respondTo ! r
      case Failure(t) => respondTo ! onFailure(t)
    }
  }

  def performCallback(): Future[Unit] = {
    val headers = List[HttpHeader]()  // TODO add auth
    for {
      entity <- Marshal(callbackMessage).to[RequestEntity]
      request = HttpRequest(method = HttpMethods.POST, uri = callbackUri.toString, entity = entity)
      // TODO add logging for retries
      response <- withRetry(() => Http().singleRequest(request.withHeaders(headers)), backoff = retryBackoff)
      result <- if (response.status.isFailure()) {
        response.entity.dataBytes.runFold(ByteString(""))(_ ++ _).map(_.utf8String) flatMap { errorBody =>
          Future.failed(
            new RuntimeException(s"Permanently failed to send callback for workflow state $lastWorkflowState to $callbackUri : $errorBody")
          )
        }
      } else {
        log.info(s"Successfully sent callback for workflow state $lastWorkflowState to $callbackUri")
        // TODO send metadata
        Future.unit
      }
    } yield result
  }
}

