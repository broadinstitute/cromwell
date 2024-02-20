package cromwell.engine.workflow.lifecycle.finalization

import akka.Done
import akka.actor.{Actor, ActorLogging, ActorRef, ActorSystem, Props}
import akka.event.LoggingReceive
import akka.http.scaladsl.Http
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.marshalling.Marshal
import akka.http.scaladsl.model._
import akka.http.scaladsl.model.headers.RawHeader
import akka.routing.Broadcast
import akka.util.ByteString
import cats.data.Validated.{Invalid, Valid}
import cats.implicits.{catsSyntaxValidatedId, toTraverseOps}
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import common.validation.ErrorOr
import cromwell.cloudsupport.azure.AzureCredentials
import cromwell.core.Dispatcher.IoDispatcher
import cromwell.core.retry.Retry.withRetry
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.core.{CallOutputs, WorkflowId, WorkflowState}
import cromwell.engine.workflow.lifecycle.finalization.WorkflowCallbackActor.PerformCallbackCommand
import cromwell.engine.workflow.lifecycle.finalization.WorkflowCallbackJsonSupport._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.{ExecutionContext, Future}
import scala.concurrent.duration.DurationInt
import scala.util.Try
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import java.net.URI
import java.time.Instant
import java.util.concurrent.Executors
import scala.util.{Failure, Success}

case class WorkflowCallbackConfig(enabled: Boolean,
                                  numThreads: Int,
                                  retryBackoff: SimpleExponentialBackoff,
                                  maxRetries: Int,
                                  defaultUri: Option[URI], // May be overridden by workflow options
                                  authMethod: Option[WorkflowCallbackConfig.AuthMethod]
)

object WorkflowCallbackConfig extends LazyLogging {
  sealed trait AuthMethod { def getAccessToken: ErrorOr.ErrorOr[String] }
  case object AzureAuth extends AuthMethod {
    override def getAccessToken: ErrorOr.ErrorOr[String] = AzureCredentials.getAccessToken()
  }

  case class StaticTokenAuth(token: String) extends AuthMethod {
    override def getAccessToken: ErrorOr.ErrorOr[String] = token.validNel
  }

  private lazy val defaultNumThreads = 5
  private lazy val defaultRetryBackoff = SimpleExponentialBackoff(3.seconds, 5.minutes, 1.1)
  private lazy val defaultMaxRetries = 10

  def empty: WorkflowCallbackConfig = WorkflowCallbackConfig(
    false,
    defaultNumThreads,
    defaultRetryBackoff,
    defaultMaxRetries,
    None,
    None
  )

  def apply(config: Config): WorkflowCallbackConfig = {
    val enabled = config.as[Boolean]("enabled")
    val numThreads = config.as[Option[Int]]("num-threads").getOrElse(defaultNumThreads)
    val backoff =
      config.as[Option[Config]]("request-backoff").map(SimpleExponentialBackoff(_)).getOrElse(defaultRetryBackoff)
    val maxRetries = config.as[Option[Int]]("max-retries").getOrElse(defaultMaxRetries)
    val uri = config.as[Option[String]]("endpoint").flatMap(createAndValidateUri)

    val authMethod = if (config.hasPath("auth.azure")) {
      Option(AzureAuth)
    } else None

    WorkflowCallbackConfig(
      enabled = enabled,
      numThreads = numThreads,
      retryBackoff = backoff,
      maxRetries = maxRetries,
      defaultUri = uri,
      authMethod = authMethod
    )
  }

  def createAndValidateUri(uriString: String): Option[URI] =
    Try(new URI(uriString)) match {
      case Success(uri) => Option(uri)
      case Failure(err) =>
        logger.warn(s"Failed to parse provided workflow callback URI (${uriString}): $err")
        None
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
                                          workflowOutputs: CallOutputs,
                                          failureMessage: List[String]
  )

  def props(serviceRegistryActor: ActorRef,
            callbackConfig: WorkflowCallbackConfig,
            httpClient: CallbackHttpHandler = CallbackHttpHandlerImpl
  ) = Props(
    new WorkflowCallbackActor(serviceRegistryActor, callbackConfig, httpClient)
  ).withDispatcher(IoDispatcher)
}

class WorkflowCallbackActor(serviceRegistryActor: ActorRef,
                            config: WorkflowCallbackConfig,
                            httpClient: CallbackHttpHandler
) extends Actor
    with ActorLogging {

  // Create a dedicated thread pool for this actor so its damage is limited if we end up with
  // too many threads all taking a long time to do callbacks. If we're frequently saturating
  // this thread pool, consider refactoring this actor to maintain a queue of callbacks and
  // handle them one at a time, as opposed to starting a thread to perform each as soon as its
  // received.
  implicit val ec = ExecutionContext.fromExecutor(Executors.newFixedThreadPool(config.numThreads))
  implicit val system: ActorSystem = context.system

  override def receive: Actor.Receive = LoggingReceive {
    case PerformCallbackCommand(workflowId, requestedCallbackUri, terminalState, outputs, failures) =>
      // If no uri was provided to us here, fall back to the one in config. If there isn't
      // one there, do not perform a callback.
      val callbackUri: Option[URI] =
        requestedCallbackUri.map(WorkflowCallbackConfig.createAndValidateUri).getOrElse(config.defaultUri)
      callbackUri
        .map { uri =>
          performCallback(workflowId, uri, terminalState, outputs, failures) onComplete {
            case Success(_) =>
              log.info(
                s"Successfully sent callback for workflow for workflow $workflowId in state $terminalState to $uri"
              )
              sendMetadata(workflowId, successful = true, uri)
            case Failure(t) =>
              log.warning(
                s"Permanently failed to send callback for workflow $workflowId in state $terminalState to $uri: ${t.getMessage}"
              )
              sendMetadata(workflowId, successful = false, uri)
          }
        }
        .getOrElse(())
    case Broadcast(ShutdownCommand) | ShutdownCommand => context stop self
    case other => log.warning(s"WorkflowCallbackActor received an unexpected message: $other")
  }

  private def makeHeaders: Future[List[HttpHeader]] =
    config.authMethod.toList
      .map(_.getAccessToken)
      .map {
        case Valid(header) => Future.successful(header)
        case Invalid(err) => Future.failed(new RuntimeException(err.toString))
      }
      .map(t => t.map(t => RawHeader("Authorization", s"Bearer $t")))
      .traverse(identity)

  private def performCallback(workflowId: WorkflowId,
                              callbackUri: URI,
                              terminalState: WorkflowState,
                              outputs: CallOutputs,
                              failures: List[String]
  ): Future[Done] = {
    val callbackPostBody = CallbackMessage(
      workflowId.toString,
      terminalState.toString,
      outputs.outputs.map(entry => (entry._1.identifier.fullyQualifiedName.value, entry._2)),
      failures
    )
    for {
      entity <- Marshal(callbackPostBody).to[RequestEntity]
      headers <- makeHeaders
      request = HttpRequest(method = HttpMethods.POST, uri = callbackUri.toString, entity = entity).withHeaders(headers)
      response <- withRetry(
        () => sendRequestOrFail(request),
        backoff = config.retryBackoff,
        maxRetries = Option(config.maxRetries),
        onRetry = err =>
          log.warning(
            s"Will retry after failure to send workflow callback for workflow $workflowId in state $terminalState to $callbackUri : $err"
          )
      )
      result <-
        // Akka will get upset if we have a response body and leave it totally unread.
        // Since there's nothing here we want to read, we need to deliberately discard it.
        response.entity.discardBytes().future
    } yield result
  }

  private def sendRequestOrFail(request: HttpRequest): Future[HttpResponse] =
    httpClient
      .sendRequest(request)
      .flatMap(response =>
        if (response.status.isFailure()) {
          response.entity.dataBytes.runFold(ByteString(""))(_ ++ _).map(_.utf8String) flatMap { errorBody =>
            Future.failed(
              new RuntimeException(s"HTTP ${response.status.value}: $errorBody")
            )
          }
        } else Future.successful(response)
      )

  private def sendMetadata(workflowId: WorkflowId, successful: Boolean, uri: URI): Unit = {
    val events = List(
      MetadataEvent(
        MetadataKey(workflowId, None, "workflowCallback", "successful"),
        MetadataValue(successful)
      ),
      MetadataEvent(
        MetadataKey(workflowId, None, "workflowCallback", "url"),
        MetadataValue(uri.toString)
      ),
      MetadataEvent(
        MetadataKey(workflowId, None, "workflowCallback", "timestamp"),
        MetadataValue(Instant.now())
      )
    )
    serviceRegistryActor ! PutMetadataAction(events)
  }

}

// Wrap the http call in a trait so it can be easily mocked for testing
trait CallbackHttpHandler {
  def sendRequest(httpRequest: HttpRequest)(implicit actorSystem: ActorSystem): Future[HttpResponse]
}

object CallbackHttpHandlerImpl extends CallbackHttpHandler {
  override def sendRequest(httpRequest: HttpRequest)(implicit actorSystem: ActorSystem): Future[HttpResponse] =
    Http().singleRequest(httpRequest)
}
