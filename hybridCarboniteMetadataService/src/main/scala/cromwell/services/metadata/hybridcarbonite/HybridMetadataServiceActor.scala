package cromwell.services.metadata.hybridcarbonite

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.event.LoggingReceive
import cats.data.NonEmptyList
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.services.metadata.MetadataService.{MetadataReadAction, MetadataServiceAction, MetadataWriteAction, ValidateWorkflowIdInMetadata}
import cromwell.services.metadata.impl.{MetadataServiceActor, ReadMetadataRegulatorActor}
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.ExecutionContext

class HybridMetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with ActorLogging with GracefulShutdownHelper {

  implicit val ec: ExecutionContext = context.dispatcher

  def getSubserviceConfig(configPath: String) = if (serviceConfig.hasPath(configPath)) serviceConfig.getConfig(configPath) else {
    log.warning(s"No ${getClass.getSimpleName} subservice config found for $configPath: assuming empty.")
    ConfigFactory.empty()
  }

  val classicMetadataServiceActorConfig = getSubserviceConfig("classicMetadataServiceActor")
  val carboniteMetadataServiceActorConfig = getSubserviceConfig("carboniteMetadataServiceActor")

  val classicMetadataService = context.actorOf(MetadataServiceActor.props(classicMetadataServiceActorConfig, globalConfig, serviceRegistryActor).withDispatcher(ServiceDispatcher))
  val carboniteMetadataService = context.actorOf(CarboniteMetadataServiceActor.props(carboniteMetadataServiceActorConfig, globalConfig, serviceRegistryActor).withDispatcher(ServiceDispatcher))

  def readDeciderActorProps(): Props = {
    HybridReadDeciderActor.props(classicMetadataService, carboniteMetadataService).withDispatcher(ServiceDispatcher)
  }

  val readRegulatorActor = context.actorOf(ReadMetadataRegulatorActor.props(readDeciderActorProps, readDeciderActorProps), "ReadMetadataRegulatorActor_for_HMSA")

  override def receive = LoggingReceive(akka.event.Logging.InfoLevel) {
    case action: MetadataServiceAction => action match {
      case read: MetadataReadAction => readRegulatorActor forward read
      case write: MetadataWriteAction => classicMetadataService.forward(write)
      case validate: ValidateWorkflowIdInMetadata => classicMetadataService forward validate
    }

    case ShutdownCommand => waitForActorsAndShutdown(NonEmptyList.of(classicMetadataService, carboniteMetadataService))
  }
}

object HybridMetadataServiceActor {
  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) = Props(new HybridMetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef))
}
