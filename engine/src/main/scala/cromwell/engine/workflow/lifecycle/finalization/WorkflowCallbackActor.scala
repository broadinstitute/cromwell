package cromwell.engine.workflow.lifecycle.finalization

import akka.actor.{Actor, ActorLogging, ActorRef, ActorSystem, Props}
import akka.event.LoggingReceive
import akka.http.scaladsl.Http
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.marshalling.Marshal
import akka.http.scaladsl.model._
import akka.http.scaladsl.model.headers.RawHeader
import akka.util.ByteString
import cats.data.Validated.{Invalid, Valid}
import cats.implicits.toTraverseOps
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import common.validation.ErrorOr
import cromwell.cloudsupport.azure.AzureCredentials
import cromwell.core.Dispatcher.IoDispatcher
import cromwell.core.retry.Retry.withRetry
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.core.{CallOutputs, WorkflowId, WorkflowState, WorkflowSucceeded}
import cromwell.engine.workflow.lifecycle.finalization.WorkflowCallbackActor.PerformCallbackCommand
import cromwell.engine.workflow.lifecycle.finalization.WorkflowCallbackConfig.AuthMethod
import cromwell.engine.workflow.lifecycle.finalization.WorkflowCallbackJsonSupport._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.ExecutionContextExecutor
import scala.concurrent.duration.DurationInt
import scala.util.Try
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}

import java.net.URI
import java.time.Instant
import scala.concurrent.Future
import scala.util.{Failure, Success}


case class WorkflowCallbackConfig(enabled: Boolean,
                                  defaultUri: Option[URI], // May be overridden by workflow options
                                  retryBackoff: Option[SimpleExponentialBackoff],
                                  maxRetries: Option[Int],
                                  authMethod: Option[WorkflowCallbackConfig.AuthMethod])

object WorkflowCallbackConfig extends LazyLogging {
  sealed trait AuthMethod { def getAccessToken:  ErrorOr.ErrorOr[String]  }
  case object AzureAuth extends AuthMethod {
    override def getAccessToken: ErrorOr.ErrorOr[String] = AzureCredentials.getAccessToken()
  }
  // TODO
  //  case class GoogleAuth(mode: GoogleAuthMode) extends AuthMethod {
  //    override def getAccessToken: String = ???
  //  }

  def apply(config: Config): WorkflowCallbackConfig = {
    val backoff = config.as[Option[Config]]("request-backoff").map(SimpleExponentialBackoff(_))
    val maxRetries = config.as[Option[Int]]("max-retries")
    val uri = config.as[Option[String]]("endpoint").flatMap(createAndValidateUri)

    // TODO add google auth
    val authMethod = if (config.hasPath("auth.azure")) {
      Option(AzureAuth)
    } else None

    WorkflowCallbackConfig(
      config.getBoolean("enabled"),
      defaultUri = uri,
      retryBackoff = backoff,
      maxRetries = maxRetries,
      authMethod = authMethod
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

  final case class PerformCallbackCommand(workflowId: WorkflowId,
                                          uri: Option[String],
                                          terminalState: WorkflowState,
                                          workflowOutputs: CallOutputs)

  def props(serviceRegistryActor: ActorRef,
            defaultCallbackUri: Option[URI],
            retryBackoff: Option[SimpleExponentialBackoff] = None,
            maxRetries: Option[Int] = None,
            authMethod: Option[WorkflowCallbackConfig.AuthMethod] = None
           ) = Props(
    new WorkflowCallbackActor(
      serviceRegistryActor,
      defaultCallbackUri,
      retryBackoff,
      maxRetries,
      authMethod)
  ).withDispatcher(IoDispatcher)
}

class WorkflowCallbackActor(serviceRegistryActor: ActorRef,
                            defaultCallbackUri: Option[URI],
                            retryBackoff: Option[SimpleExponentialBackoff],
                            maxRetries: Option[Int],
                            authMethod: Option[AuthMethod])
  extends Actor with ActorLogging {

  implicit val ec: ExecutionContextExecutor = context.dispatcher
  implicit val system: ActorSystem = context.system

  private lazy val defaultRetryBackoff = SimpleExponentialBackoff(3.seconds, 5.minutes, 1.1)
  private lazy val defaultMaxRetries = 10

  private lazy val backoffWithDefault = retryBackoff.getOrElse(defaultRetryBackoff)
  private lazy val maxRetriesWithDefault = maxRetries.getOrElse(defaultMaxRetries)

  override def receive: Actor.Receive = LoggingReceive {
    case PerformCallbackCommand(workflowId, uri, terminalState, outputs) =>
      // If no uri was provided to us here, fall back to the one in config. If there isn't
      // one there, do not perform a callback.
      val callbackUri: Option[URI] = uri.map(WorkflowCallbackConfig.createAndValidateUri).getOrElse(defaultCallbackUri)
      callbackUri.map { uri =>
        performCallback(workflowId, uri, terminalState, outputs) onComplete {
          case Success(_) =>
            log.info(s"Successfully sent callback for workflow for workflow $workflowId in state $terminalState to $uri")
            sendMetadata(workflowId, successful = true)
          case Failure(t) =>
            log.warning(s"Permanently failed to send callback for workflow $workflowId in state $terminalState to $uri: ${t.getMessage}")
            sendMetadata(workflowId, successful = false)
        }
      }.getOrElse(())
    case other => log.warning(s"WorkflowCallbackActor received an unexpected message: $other")
  }

  private def makeHeaders: Future[List[HttpHeader]] = {
    authMethod.toList.map(_.getAccessToken).map {
      case Valid(header) => Future.successful(header)
      case Invalid(err) => Future.failed(new RuntimeException(err.toString)) // TODO better error
    }
      .map(t => t.map(RawHeader("Authorization", _)))
      .traverse(identity)
  }

  private def performCallback(workflowId: WorkflowId, callbackUri: URI, terminalState: WorkflowState, outputs: CallOutputs): Future[Unit] = {
    // Only send outputs if the workflow succeeded
    val callbackPostBody = terminalState match {
      case WorkflowSucceeded => CallbackMessage(workflowId.toString, terminalState.toString, Option(outputs.outputs.map(entry => (entry._1.name, entry._2))))
      case _ => CallbackMessage(workflowId.toString, terminalState.toString, None)
    }
    for {
      entity <- Marshal(callbackPostBody).to[RequestEntity]
      headers <- makeHeaders
      request = HttpRequest(method = HttpMethods.POST, uri = callbackUri.toString, entity = entity).withHeaders(headers)
      response <- withRetry(
        () => sendRequestOrFail(request),
        backoff = backoffWithDefault,
        maxRetries = Option(maxRetriesWithDefault),
        onRetry = err => log.warning(s"Will retry after failure to send workflow callback for workflow $workflowId in state $terminalState to $defaultCallbackUri : $err")
      )
      result <- {
        response.entity.discardBytes()
        Future.unit
      }
    } yield result
  }

  private def sendRequestOrFail(request: HttpRequest): Future[HttpResponse] =
    Http().singleRequest(request).flatMap(response =>
      if (response.status.isFailure()) {
        response.entity.dataBytes.runFold(ByteString(""))(_ ++ _).map(_.utf8String) flatMap { errorBody =>
          Future.failed(
            new RuntimeException(s"HTTP ${response.status.value}: $errorBody")
          )
        }
      } else Future.successful(response)
    )

  private def sendMetadata(workflowId: WorkflowId, successful: Boolean): Unit = {
    val events = List(
      MetadataEvent(
        MetadataKey(workflowId, None, "workflowCallback", "successful"),
        MetadataValue(successful)
      ),
      MetadataEvent(
        MetadataKey(workflowId, None, "workflowCallback", "url"),
        MetadataValue(defaultCallbackUri.toString)
      ),
      MetadataEvent(
        MetadataKey(workflowId, None, "workflowCallback", "timestamp"),
        MetadataValue(Instant.now())
      )
    )
    serviceRegistryActor ! PutMetadataAction(events)
  }

}

