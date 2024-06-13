package cromwell.services.metrics.bard.impl

import akka.actor.{Actor, ActorRef, ActorSystem, Props}
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import common.util.StringUtil.EnhancedToStringable
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.services.metrics.bard.BardEventing.BardEventRequest
import cromwell.services.metrics.bard.{BardConfig, BardService}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.ExecutionContext

class BardEventingActor(serviceConfig: Config, globalConfig: Config, serviceRegistry: ActorRef)
    extends Actor
    with LazyLogging {

  implicit val system: ActorSystem = context.system
  implicit val ec: ExecutionContext = context.dispatcher

  lazy val bardConfig: BardConfig = BardConfig(serviceConfig)
  lazy val bardService: BardService = new BardService(bardConfig.baseUrl, bardConfig.connectionPoolSize, serviceRegistry)

  override def receive: Receive = {
    case BardEventRequest(event) if bardConfig.enabled => bardService.sendEvent(event)
    // This service currently doesn't do any work on shutdown but the service registry pattern requires it
    // (see https://github.com/broadinstitute/cromwell/issues/2575)
    case ShutdownCommand => context stop self
    case other =>
      logger.error(
        s"Programmer Error: Unexpected message ${other.toPrettyElidedString(1000)} received by ${this.self.path.name}."
      )
  }
}

object BardEventingActor {

  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef): Props =
    Props(new BardEventingActor(serviceConfig, globalConfig, serviceRegistryActor))
      .withDispatcher(ServiceDispatcher)
      .withMailbox("akka.bard-actor-mailbox")

}
