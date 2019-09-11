package cromwell.services.metadata.impl.carboniter

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import com.typesafe.config.Config
import cromwell.services.metadata.MetadataService.{MetadataReadAction, MetadataWriteAction, MetadataWriteFailure}
import net.ceedubs.ficus.Ficus._
import scala.concurrent.duration._

class CarboniteMetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with ActorLogging {

  val carboniteWorker: Option[ActorRef] = {

    val carboniteInterval = serviceConfig.getOrElse[Duration]("metadata-summary-refresh-interval", default = Duration.Inf)
    if (carboniteInterval.isFinite()) {
      val finiteCarboniteInterval = carboniteInterval.asInstanceOf[FiniteDuration]
      Option(context.actorOf(CarboniteWorkerActor.props(finiteCarboniteInterval)))
    } else {
      log.info("Carboniting interval is not specified. No metadata carboniting will be performed.")
      None
    }
  }

  override def receive: Receive = {
    case read: MetadataReadAction =>
      val worker = context.actorOf(CarbonitedMetadataReaderActor.props)
      worker forward read
    case write: MetadataWriteAction =>
      val error = new NotImplementedError(s"Programmer Error! Carboniter Worker should never be sent write requests (but got $write from $sender)")
      sender ! MetadataWriteFailure(error, write.events)
  }
}

object CarboniteMetadataServiceActor {
  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) = Props(new CarboniteMetadataServiceActor(serviceConfig, globalConfig, serviceRegistryActor))
}
