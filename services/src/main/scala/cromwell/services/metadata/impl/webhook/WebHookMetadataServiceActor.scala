package cromwell.services.metadata.impl.webhook

import akka.actor.{Actor, ActorLogging, ActorRef, ActorSystem, Props}
import akka.http.scaladsl.Http
import akka.http.scaladsl.model._
import com.typesafe.config.Config
import cromwell.core.Dispatcher._
import cromwell.services.metadata.MetadataService.{MetadataWriteFailure, MetadataWriteSuccess, PutMetadataAction, PutMetadataActionAndRespond}
import cromwell.services.metadata._
import cromwell.services.metadata.impl.MetadataServiceActor
import net.ceedubs.ficus.Ficus._
import spray.json._

import scala.concurrent.{ExecutionContextExecutor, Future}
import scala.util.{Failure, Success}


class WebHookHybridMetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with ActorLogging {

  val standardMetadataActor: ActorRef = context.actorOf(MetadataServiceActor.props(serviceConfig, globalConfig, serviceRegistryActor))
  val webHookActor: ActorRef = context.actorOf(WebHookMetadataServiceActor.props(serviceConfig, globalConfig, serviceRegistryActor))

  override def receive = {
    case action: PutMetadataAction =>
      standardMetadataActor forward action
      webHookActor forward action
    case action: PutMetadataActionAndRespond =>
      standardMetadataActor forward action
      webHookActor forward PutMetadataAction(action.events)
    case anythingElse => standardMetadataActor forward anythingElse
  }
}

/**
  * A *write-only* metadata service implementation which pushes all events to an HTTP webhook.
  * The expectation is that metadata reads are being handled outside of this Cromwell instance.
  */
class WebHookMetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with ActorLogging {
  implicit val ec: ExecutionContextExecutor = context.dispatcher
  implicit val system: ActorSystem = context.system
//  implicit val materializer = ActorMaterializer(system)

  val httphook: String = serviceConfig.as[String]("url")

  override def receive = {
    case action: PutMetadataAction =>
      publishMessages(action.events).failed foreach { e =>
        log.error(e, "Failed to post metadata: " + action.events)
      }
    case action: PutMetadataActionAndRespond =>
      publishMessages(action.events) onComplete {
        case Success(_) => action.replyTo ! MetadataWriteSuccess(action.events)
        case Failure(e) => action.replyTo ! MetadataWriteFailure(e, action.events)
      }
  }


  private def publishMessages(events: Iterable[MetadataEvent]): Future[Unit] = {

    import WebHookMetadataServiceActor.EnhancedMetadataEvents

    val eventsJson = events.toJson
    log.debug("Publishing metadata event to " + httphook)

    // Todo, turn this into a proper JSON object
    val entity = HttpEntity(ContentTypes.`application/json`, eventsJson.mkString(","))

    val httpRequest = HttpRequest(
      method = HttpMethods.POST,
      uri = httphook,
      headers = List[HttpHeader](),
      entity=entity
    )


    val responseFuture = Http().singleRequest(httpRequest)

    Future {
      responseFuture
          .onComplete {
            case Success(_) => None
            case Failure(exception) => log.error("Failed to send webhook metadata: " + exception.toString)
          }
    }
  }
}

object WebHookMetadataServiceActor {
  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) = {
    Props(new WebHookMetadataServiceActor(serviceConfig, globalConfig, serviceRegistryActor)).withDispatcher(ServiceDispatcher)
  }

  implicit class EnhancedMetadataEvents(val e: Iterable[MetadataEvent]) extends AnyVal {
    import MetadataJsonSupport._

    def toJson: Seq[String] = e.map(_.toJson.toString()).toSeq
  }
}

