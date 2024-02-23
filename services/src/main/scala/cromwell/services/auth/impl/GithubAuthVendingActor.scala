package cromwell.services.auth.impl

import akka.actor.{Actor, ActorRef, Props}
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.services.auth.GithubAuthVending.{GithubAuthRequest, GithubAuthTokenResponse, NoGithubAuthResponse}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

class GithubAuthVendingActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef)
    extends Actor
    with LazyLogging {

  lazy val enabled = serviceConfig.getBoolean("enabled")

  override def receive: Receive = {
    case GithubAuthRequest(_) if enabled =>
      sender() ! GithubAuthTokenResponse(serviceConfig.getString("access-token"))
    // This service currently doesn't do any work on shutdown but the service registry pattern requires it
    // (see https://github.com/broadinstitute/cromwell/issues/2575)
    case ShutdownCommand => context stop self
    case _ =>
      sender() ! NoGithubAuthResponse
  }
}

object GithubAuthVendingActor {

  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef): Props =
    Props(new GithubAuthVendingActor(serviceConfig, globalConfig, serviceRegistryActor))
      .withDispatcher(ServiceDispatcher)

}
