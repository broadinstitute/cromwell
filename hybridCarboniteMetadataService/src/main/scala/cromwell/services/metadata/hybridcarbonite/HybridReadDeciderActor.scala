package cromwell.services.metadata.hybridcarbonite

import akka.actor.{Actor, ActorRef, PoisonPill, Props}
import akka.event.LoggingReceive
import cromwell.services.metadata.MetadataService.{MetadataReadAction, MetadataServiceResponse}
import cromwell.services.metadata.impl.builder.MetadataBuilderActor.FailedMetadataResponse

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

class HybridReadDeciderActor(classicMetadataServiceActor: ActorRef, carboniteMetadataServiceActor: ActorRef) extends Actor {

  implicit val ec: ExecutionContext = context.dispatcher

  // TODO: Better decision logic!
  def decide(read: MetadataReadAction): Future[ActorRef] = Future.successful(classicMetadataServiceActor)

  var sndr: Option[ActorRef] = None

  override def receive: Receive = LoggingReceive(akka.event.Logging.InfoLevel) {
    case read: MetadataReadAction =>
      sndr = Option(sender())
      decide(read) onComplete {
        case Success(forwardTo) =>
          forwardTo ! read
        case Failure(e) =>
          respondAndStop(FailedMetadataResponse(read, new Exception("Failed to decide which metadata service to forward request to", e)))
      }
    case response: MetadataServiceResponse =>
      respondAndStop(response)
  }

  def respondAndStop(msg: Any) = {
    sndr.foreach { _ ! msg }
    self ! PoisonPill
  }
}

object HybridReadDeciderActor {
  def props(classicMetadataServiceActor: ActorRef, carboniteMetadataServiceActor: ActorRef) =
    Props(new HybridReadDeciderActor(classicMetadataServiceActor, carboniteMetadataServiceActor))
}
