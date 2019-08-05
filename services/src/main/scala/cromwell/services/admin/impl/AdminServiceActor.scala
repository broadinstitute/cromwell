package cromwell.services.admin.impl

import akka.actor.{Actor, ActorRef, Props}
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.services.admin.AdminServiceMessages.{ListSubmissionsRequest, PauseSubmissionRequest, PauseSubmissionSuccess}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.ExecutionContext

class AdminServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with LazyLogging {

  implicit val ec: ExecutionContext = context.dispatcher

  override def receive: Receive = {
    case PauseSubmissionRequest => sender ! PauseSubmissionSuccess
    case ListSubmissionsRequest => sender ! PauseSubmissionSuccess

    case ShutdownCommand =>
      // This service doesn't do any work on shutdown but the service registry pattern requires it (see #2575)
      context.stop(self)
  }

}

object AdminServiceActor {
  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) =
    Props(new AdminServiceActor(serviceConfig, globalConfig, serviceRegistryActor)).withDispatcher(ServiceDispatcher)
}
