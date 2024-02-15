package cromwell.services.auth.impl

import akka.actor.{Actor, ActorRef, Props}
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.services.auth.GithubAuthVending.{GithubAuthRequest, GithubAuthTokenResponse, NoGithubAuthResponse}

class GithubAuthVendingActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef)
    extends Actor
    with LazyLogging {

  lazy val enabled = serviceConfig.getBoolean("enabled")

  override def receive: Receive = {
    case GithubAuthRequest(_) if enabled =>
      sender() ! GithubAuthTokenResponse(serviceConfig.getString("access-token"))
    case _ =>
      sender() ! NoGithubAuthResponse
  }
}

object GithubAuthVendingActor {

  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef): Props =
    Props(new GithubAuthVendingActor(serviceConfig, globalConfig, serviceRegistryActor))
      .withDispatcher(ServiceDispatcher)

}
