package cromwell.services.metadata.hybridcarbonite

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.services.metadata.MetadataService.{MetadataReadAction, MetadataServiceAction, MetadataWriteAction, ValidateWorkflowIdInMetadata}
import cromwell.services.metadata.impl.MetadataServiceActor
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

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

  // Determines which sub-actor to forward the read query to
  def whereToForward(metadataReadAction: MetadataReadAction): Future[ActorRef] = {
    // For now, this is hard coded to the classic metadata service:
    Future.successful(classicMetadataService)
  }

  override def receive = {
    case action: MetadataServiceAction => action match {
      case read: MetadataReadAction =>
        val sndr = sender()
        whereToForward(read) onComplete {
          case Success(metadataService: ActorRef) => metadataService.tell(read, sndr)
          case Failure(exception) =>
            // Log the error, and fall back to the classic actor:
            log.error(exception, s"Failed to decide which metadata service sub-actor to send $read to. Falling back to classic")
            classicMetadataService.tell(read, sndr)
        }

      case write: MetadataWriteAction => classicMetadataService.forward(write)
      case validate: ValidateWorkflowIdInMetadata => classicMetadataService forward validate
    }



    case ShutdownCommand => waitForActorsAndShutdown(NonEmptyList.of(classicMetadataService, carboniteMetadataService))
  }
}

object HybridMetadataServiceActor {
  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) = Props(new HybridMetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef))
}
