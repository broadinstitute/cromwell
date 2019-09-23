package cromwell.services.metadata.hybridcarbonite

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import cromwell.services.metadata.MetadataService.{MetadataReadAction, MetadataWriteAction, MetadataWriteFailure}
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

class CarboniteMetadataServiceActor(carboniteConfig: HybridCarboniteConfig, serviceRegistryActor: ActorRef) extends Actor with ActorLogging with GracefulShutdownHelper {

  // TODO: [CARBONITE] Pass in a real reference to the IO actor somehow
  val ioActor = context.actorOf(Props(new Actor {
    override def receive: Receive = Actor.emptyBehavior }))

  // TODO: [CARBONITE] Validate the carbonite config
  def thawingActorProps(): Props = CarbonitedMetadataThawingActor.props(carboniteConfig, serviceRegistryActor, ioActor)

  val carboniteWorker: Option[ActorRef] = {
    carboniteConfig.carboniteInterval map { interval => context.actorOf(CarboniteWorkerActor.props(interval)) }
  }

  override def receive: Receive = {
    case read: MetadataReadAction =>
      val worker = context.actorOf(CarbonitedMetadataReaderActor.props)
      worker forward read
    case write: MetadataWriteAction =>
      val error = new UnsupportedOperationException(s"Programmer Error! Carboniter Worker should never be sent write requests (but got $write from $sender)")
      sender ! MetadataWriteFailure(error, write.events)
    case ShutdownCommand =>
      carboniteWorker match {
        case Some(worker) => waitForActorsAndShutdown(NonEmptyList.of(worker))
        case None => context.stop(self)
      }
  }
}

object CarboniteMetadataServiceActor {
  def props(serviceConfig: HybridCarboniteConfig, serviceRegistryActor: ActorRef) = Props(new CarboniteMetadataServiceActor(serviceConfig, serviceRegistryActor))
}
