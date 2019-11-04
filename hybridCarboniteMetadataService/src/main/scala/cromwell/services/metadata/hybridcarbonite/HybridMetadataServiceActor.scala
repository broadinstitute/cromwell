package cromwell.services.metadata.hybridcarbonite

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import com.typesafe.config.Config
import common.exception.AggregatedMessageException
import common.validation.Checked._
import common.validation.Validation._
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.services.metadata.MetadataService.{BuildMetadataJsonAction, MetadataServiceAction, MetadataWriteAction, ValidateWorkflowIdInMetadata}
import cromwell.services.metadata.impl.{MetadataServiceActor, ReadMetadataRegulatorActor}
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import cromwell.services.metadata.hybridcarbonite.HybridMetadataServiceActor.CarboniteConfigPath

import scala.concurrent.ExecutionContext

class HybridMetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with ActorLogging with GracefulShutdownHelper {

  implicit val ec: ExecutionContext = context.dispatcher

  val carboniteMetadataServiceActorConfig: HybridCarboniteConfig = {
    val hybridCarboniteConfig = if (serviceConfig.hasPath(CarboniteConfigPath)) {
      HybridCarboniteConfig.parseConfig(serviceConfig.getConfig(CarboniteConfigPath))(context.system).contextualizeErrors("parse config")
    } else {
      s"No ${getClass.getSimpleName} subservice config found for $CarboniteConfigPath.".invalidNelCheck
    }

    hybridCarboniteConfig match {
      case Right(config) => config
      case Left(e) => throw AggregatedMessageException("Failed to initialize HybridMetadataServiceActor", e.toList)
    }
  }

  val classicMetadataService = context.actorOf(MetadataServiceActor.props(serviceConfig, globalConfig, serviceRegistryActor).withDispatcher(ServiceDispatcher))
  val carboniteMetadataService = context.actorOf(CarboniteMetadataServiceActor.props(carboniteMetadataServiceActorConfig, serviceRegistryActor).withDispatcher(ServiceDispatcher))

  def readDeciderActorProps(): Props = {
    HybridReadDeciderActor.props(classicMetadataService, carboniteMetadataService).withDispatcher(ServiceDispatcher)
  }

  val readRegulatorActor = context.actorOf(ReadMetadataRegulatorActor.props(readDeciderActorProps, readDeciderActorProps), "ReadMetadataRegulatorActor_for_HMSA")

  override def receive = {
    case action: MetadataServiceAction => action match {
      case read: BuildMetadataJsonAction => readRegulatorActor forward read
      case write: MetadataWriteAction => classicMetadataService.forward(write)
      case validate: ValidateWorkflowIdInMetadata => classicMetadataService forward validate
    }
    case ShutdownCommand => waitForActorsAndShutdown(NonEmptyList.of(classicMetadataService, carboniteMetadataService))
  }
}

object HybridMetadataServiceActor {
  val CarboniteConfigPath = "carbonite-metadata-service"

  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) = Props(new HybridMetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef))
}
