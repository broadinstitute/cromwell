package cromwell.services.metadata.impl

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import com.typesafe.config.Config
import cromwell.services.metadata.MetadataService.{PutMetadataAction, PutMetadataActionAndRespond}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import org.slf4j.LoggerFactory
import cromwell.core.Dispatcher._

class LogSplittingMetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with ActorLogging {

  val logger = LoggerFactory.getLogger("metadata_logger")
  val logMetadata: ActorRef = context.actorOf(Props(LogMetadataActor(logger.info(_))).withDispatcher(ServiceDispatcher))

  val standardMetadataActor: ActorRef = context.actorOf(MetadataServiceActor.props(serviceConfig, globalConfig, serviceRegistryActor))

  override def receive = {
    case action: PutMetadataAction =>
      standardMetadataActor forward action
      logMetadata forward action
    case action: PutMetadataActionAndRespond =>
      standardMetadataActor forward action
      logMetadata forward PutMetadataAction(action.events)
    case sc@ShutdownCommand =>
      context stop logMetadata
      standardMetadataActor forward sc
      context stop self
    case anythingElse => standardMetadataActor forward anythingElse
  }
}
