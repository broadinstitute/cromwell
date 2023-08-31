package cromwell.engine.workflow.lifecycle.finalization

import akka.actor.{Actor, ActorLogging, ActorRef, ActorSystem, Props}
import akka.event.LoggingReceive
import akka.http.scaladsl.Http
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.marshalling.Marshal
import akka.http.scaladsl.model._
import akka.util.ByteString
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import cromwell.backend.BackendLifecycleActor.BackendWorkflowLifecycleActorResponse
import cromwell.backend.BackendWorkflowFinalizationActor.{FinalizationSuccess, Finalize}
import cromwell.core.Dispatcher.IoDispatcher
import cromwell.core.retry.Retry.withRetry
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.core.{CallOutputs, WorkflowId, WorkflowState}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.finalization.WorkflowCallbackJsonSupport._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.ExecutionContextExecutor
import scala.concurrent.duration.DurationInt
import scala.util.Try
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}

import java.net.URI
import java.time.Instant
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}


case class WorkflowCallbackConfig(enabled: Boolean,
                                  uri: Option[URI], // May be overridden by workflow options
                                  retryBackoff: Option[SimpleExponentialBackoff]
                                  // TODO auth
  )

object WorkflowCallbackConfig extends LazyLogging {
  def apply(config: Config): WorkflowCallbackConfig = {
    val backoff = config.as[Option[Config]]("request-backoff").map(SimpleExponentialBackoff(_))
    val uri = config.as[Option[String]]("endpoint").flatMap(createAndValidateUri)

    WorkflowCallbackConfig(
      config.getBoolean("enabled"),
      uri = uri,
      retryBackoff = backoff
    )
  }

  def createAndValidateUri(uriString: String): Option[URI] = {
    // TODO validate
    Try(new URI(uriString)) match {
      case Success(uri) => Option(uri)
      case Failure(err) =>
        logger.warn(s"Failed to parse provided workflow callback URI (${uriString}): $err")
        None
    }
  }
}


/**
  * The WorkflowCallbackActor is responsible for sending a message on workflow completion to a configured endpoint.
  * This allows for users to build automation around workflow completion without needing to poll Cromwell for
  * workflow status.
  */
object WorkflowCallbackActor {
  def props(workflowId: WorkflowId,
            serviceRegistryActor: ActorRef,
            workflowDescriptor: EngineWorkflowDescriptor,
            workflowOutputs: CallOutputs,
            lastWorkflowState: WorkflowState,
            callbackUri: URI,
            retryBackoff: Option[SimpleExponentialBackoff] = None
           ) = Props(
    new WorkflowCallbackActor(
      workflowId,
      serviceRegistryActor,
      workflowDescriptor,
      workflowOutputs,
      lastWorkflowState,
      callbackUri,
      retryBackoff)
  ).withDispatcher(IoDispatcher)
}

class WorkflowCallbackActor(workflowId: WorkflowId,
                            serviceRegistryActor: ActorRef,
                            val workflowDescriptor: EngineWorkflowDescriptor,
                            workflowOutputs: CallOutputs,
                            lastWorkflowState: WorkflowState,
                            callbackUri: URI,
                            retryBackoff: Option[SimpleExponentialBackoff])
  extends Actor with ActorLogging {

  implicit val ec: ExecutionContextExecutor = context.dispatcher
  implicit val system: ActorSystem = context.system

  lazy private val callbackMessage = CallbackMessage(
    workflowId.toString, lastWorkflowState.toString, workflowOutputs.outputs.map(entry => (entry._1.name, entry._2))
  )

  private lazy val defaultRetryBackoff = SimpleExponentialBackoff(3.seconds, 5.minutes, 1.1)
  private lazy val backoffWithDefault = retryBackoff.getOrElse(defaultRetryBackoff)

  override def receive: Actor.Receive = LoggingReceive {
    case Finalize => performActionThenRespond(
      performCallback().map( _ => FinalizationSuccess),
      _ => FinalizationSuccess // Don't fail the workflow when this step fails
    )
  }

  private def performActionThenRespond(operation: => Future[BackendWorkflowLifecycleActorResponse],
                                       onFailure: (Throwable) => BackendWorkflowLifecycleActorResponse)
                                      (implicit ec: ExecutionContext): Unit = {
    val respondTo: ActorRef = sender()
    operation onComplete {
      case Success(r) =>
        log.info(s"Successfully sent callback for workflow for workflow $workflowId in state $lastWorkflowState to $callbackUri")
        sendMetadata(successful = true)
        respondTo ! r
      case Failure(t) =>
        log.warning(s"Permanently failed to send callback for workflow $workflowId in state $lastWorkflowState to $callbackUri: ${t.getMessage}")
        sendMetadata(successful = false)
        respondTo ! onFailure(t)
    }
  }

  def performCallback(): Future[Unit] = {
    val headers = List[HttpHeader]()  // TODO add auth
    for {
      entity <- Marshal(callbackMessage).to[RequestEntity]
      request = HttpRequest(method = HttpMethods.POST, uri = callbackUri.toString, entity = entity)
      response <- withRetry(
        () => sendRequestOrFail(request.withHeaders(headers)),
        backoff = backoffWithDefault,
        onRetry = err => log.warning(s"Will retry after failure to send workflow callback for workflow $workflowId in state $lastWorkflowState to $callbackUri : $err")
      )
      result <- {
        response.entity.discardBytes()
        Future.unit
      }
    } yield result
  }

  def sendRequestOrFail(request: HttpRequest): Future[HttpResponse] =
    Http().singleRequest(request).flatMap(response =>
      if (response.status.isFailure()) {
        response.entity.dataBytes.runFold(ByteString(""))(_ ++ _).map(_.utf8String) flatMap { errorBody =>
          Future.failed(
            new RuntimeException(errorBody)
          )
        }
      } else Future.successful(response)
    )

  def sendMetadata(successful: Boolean): Unit = {
    val events = List(
      MetadataEvent(
        MetadataKey(workflowId, None, "workflowCallback", "successful"),
        MetadataValue(successful)
      ),
      MetadataEvent(
        MetadataKey(workflowId, None, "workflowCallback", "url"),
        MetadataValue(callbackUri.toString)
      ),
      MetadataEvent(
        MetadataKey(workflowId, None, "workflowCallback", "timestamp"),
        MetadataValue(Instant.now())
      )
    )
    serviceRegistryActor ! PutMetadataAction(events)
  }

}

