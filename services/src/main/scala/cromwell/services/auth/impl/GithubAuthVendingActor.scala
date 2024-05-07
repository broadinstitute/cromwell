package cromwell.services.auth.impl

import akka.actor.{Actor, ActorRef, ActorSystem, Props}
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import common.util.StringUtil.EnhancedToStringable
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.services.auth.GithubAuthVending.{
  GithubAuthRequest,
  GithubAuthTokenResponse,
  GithubAuthVendingFailure,
  NoGithubAuthResponse
}
import cromwell.services.auth.ecm.{EcmConfig, EcmService}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}

class GithubAuthVendingActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef)
    extends Actor
    with LazyLogging {

  implicit val system: ActorSystem = context.system
  implicit val ec: ExecutionContext = context.dispatcher

  lazy val enabled: Boolean = serviceConfig.getBoolean("enabled")

  lazy val ecmConfigOpt: EcmConfig = EcmConfig(serviceConfig)
  lazy val ecmServiceOpt: Option[EcmService] = ecmConfigOpt.baseUrl.map(url => new EcmService(url))

  override def receive: Receive = {
    case GithubAuthRequest(terraToken) if enabled =>
      val respondTo = sender()
      ecmServiceOpt match {
        case Some(ecmService) =>
          ecmService.getGithubAccessToken(terraToken) onComplete {
            case Success(githubToken) => respondTo ! GithubAuthTokenResponse(githubToken)
            case Failure(e) => respondTo ! GithubAuthVendingFailure(e.getMessage)
          }
        case None =>
          respondTo ! GithubAuthVendingFailure(
            "Invalid configuration for service 'GithubAuthVending': missing 'ecm.base-url' value."
          )
      }
    case GithubAuthRequest(_) => sender() ! NoGithubAuthResponse
    // This service currently doesn't do any work on shutdown but the service registry pattern requires it
    // (see https://github.com/broadinstitute/cromwell/issues/2575)
    case ShutdownCommand => context stop self
    case other =>
      logger.error(
        s"Programmer Error: Unexpected message ${other.toPrettyElidedString(1000)} received by ${this.self.path.name}."
      )
      sender() ! GithubAuthVendingFailure(s"Received unexpected message ${other.toPrettyElidedString(1000)}.")
  }
}

object GithubAuthVendingActor {

  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef): Props =
    Props(new GithubAuthVendingActor(serviceConfig, globalConfig, serviceRegistryActor))
      .withDispatcher(ServiceDispatcher)

}
