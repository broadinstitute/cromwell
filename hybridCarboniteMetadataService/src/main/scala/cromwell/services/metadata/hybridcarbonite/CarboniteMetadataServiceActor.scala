package cromwell.services.metadata.hybridcarbonite

import java.util.concurrent.atomic.AtomicLong

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import cromwell.services.FailedMetadataJsonResponse
import cromwell.services.ServiceRegistryActor.{IoActorRef, NoIoActorRefAvailable, RequestIoActorRef}
import cromwell.services.metadata.MetadataService.{BuildWorkflowMetadataJsonAction, MetadataWriteAction, MetadataWriteFailure}
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import mouse.boolean._

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

class CarboniteMetadataServiceActor(carboniteConfig: HybridCarboniteConfig, serviceRegistryActor: ActorRef) extends Actor with ActorLogging with GracefulShutdownHelper {

  implicit val ec: ExecutionContext = context.dispatcher

  var ioActorOption: Option[ActorRef] = None
  scheduleIoActorLookup()

  def thawingActorProps(ioActor: ActorRef): Props = CarbonitedMetadataThawingActor.props(carboniteConfig, serviceRegistryActor, ioActor)

  var carboniteWorker: Option[ActorRef] = None

  override def receive: Receive = {
    case read: BuildWorkflowMetadataJsonAction =>
      ioActorOption match {
        case Some(ioActor) =>
          val worker = context.actorOf(thawingActorProps(ioActor), s"ThawingActor-${saltNumber.incrementAndGet()}-for-${read.workflowId}")
          worker forward read
        case None =>
          sender ! FailedMetadataJsonResponse(read, new Exception("Cannot create CarbonitedMetadataThawingActor: no IoActor reference available"))
      }
    case write: MetadataWriteAction =>
      val error = new UnsupportedOperationException(s"Programmer Error! Carboniter Worker should never be sent write requests (but got $write from $sender)")
      sender ! MetadataWriteFailure(error, write.events)
    case IoActorRef(ref) =>
      log.info(s"${getClass.getSimpleName} has received an IoActor reference")
      ioActorOption = Option(ref)
      carboniteWorker = carboniteConfig.enabled.option(context.actorOf(CarboniteWorkerActor.props(carboniteConfig, serviceRegistryActor, ref)))
    case NoIoActorRefAvailable =>
      log.warning(s"${getClass.getSimpleName} is still waiting for an IoActor reference")
      scheduleIoActorLookup()
    case ShutdownCommand =>
      carboniteWorker match {
        case Some(worker) => waitForActorsAndShutdown(NonEmptyList.of(worker))
        case None => context.stop(self)
      }
  }

  private def scheduleIoActorLookup(): Unit = {
    context.system.scheduler.scheduleOnce(1.second) {
      serviceRegistryActor ! RequestIoActorRef
    }
    ()
  }

  // For naming:
  // - Makes sure we don't accidentally duplicate sub-actor names
  // - We don't actually need the thread-safety, this class just has a convenient incrementAndGet method
  val saltNumber = new AtomicLong(0)
}

object CarboniteMetadataServiceActor {
  def props(serviceConfig: HybridCarboniteConfig, serviceRegistryActor: ActorRef) = Props(new CarboniteMetadataServiceActor(serviceConfig, serviceRegistryActor))
}
