package cromwell.services.auth.impl

import akka.actor.{Actor, ActorRef, Props}
import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.services.auth.GithubAuthVending.{
  GithubAuthRequest,
  GithubAuthTokenResponse,
  GithubAuthVendingFailure,
  NoGithubAuthResponse
}
import cromwell.services.auth.externalcreds.{EcmConfig, HttpEcmApiClientProvider}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

class GithubAuthVendingActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef)
    extends Actor
    with LazyLogging {

  lazy val enabled = serviceConfig.getBoolean("enabled")

  private lazy val ecmConfigOpt: Option[EcmConfig] = EcmConfig.apply(serviceConfig)
  private lazy val ecmClientProviderOpt: Option[HttpEcmApiClientProvider] =
    ecmConfigOpt.map(ecmConfig => new HttpEcmApiClientProvider(ecmConfig.baseUrl))

  override def receive: Receive = {
    case GithubAuthRequest(userToken) if enabled =>
      ecmClientProviderOpt match {
        case Some(ecmClientProvider) =>
          val ecmOauthApi = ecmClientProvider.getOauthApi(userToken)
          ecmOauthApi.getGithubAccessToken match {
            case Valid(githubToken) => sender() ! GithubAuthTokenResponse(githubToken)
            case Invalid(e) => sender() ! GithubAuthVendingFailure(e.head)
          }
        case None =>
          sender() ! GithubAuthVendingFailure(
            "Invalid configuration for service 'GithubAuthVending': missing 'ecm.base-url' value."
          )
      }
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
