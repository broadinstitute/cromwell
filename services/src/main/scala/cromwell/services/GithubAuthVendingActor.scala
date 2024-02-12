package cromwell.services

import akka.actor.{Actor, ActorRef, Props}
import akka.pattern.AskSupport
import akka.util.Timeout
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.languages.util.ImportResolver.ImportAuthProvider
import cromwell.services.GithubAuthVendingActor.GithubAuthRequest
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage

import scala.concurrent.{ExecutionContext, Future}
import scala.concurrent.duration.{Duration, DurationInt}

class GithubAuthVendingActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with LazyLogging {

  override def receive: Receive = {
    case GithubAuthRequest(terraToken, replyTo) =>
      replyTo ! GithubAuthVendingActor.GithubAuthVendingSuccess("access-token")
  }
}

object GithubAuthVendingActor {

  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef): Props =
    Props(new GithubAuthVendingActor(serviceConfig, globalConfig, serviceRegistryActor))
      .withDispatcher(ServiceDispatcher)

  sealed trait GithubAuthVendingMessage extends ServiceRegistryMessage {
    override def serviceName: String = "GithubAuthVending"
  }

  case class GithubAuthRequest(terraToken: String, replyTo: ActorRef) extends GithubAuthVendingMessage

  sealed trait GithubAuthVendingResponse extends GithubAuthVendingMessage
  case class GithubAuthVendingSuccess(accessToken: String) extends GithubAuthVendingResponse
  case class GithubAuthVendingFailure(error: Exception) extends GithubAuthVendingResponse

  trait GithubAuthVendingSupport extends AskSupport {
    def serviceRegistry: ActorRef
    implicit val timeout: Timeout = 10.seconds
    implicit val ec: ExecutionContext

    def importAuthProvider(token: String): ImportAuthProvider = new ImportAuthProvider {
      override def validHosts: List[String] = List("github.com")
      override def authHeader(): Future[Map[String, String]] = {
        serviceRegistry.ask(replyTo => GithubAuthRequest(token, replyTo)).mapTo[GithubAuthVendingResponse].flatMap {
          case GithubAuthVendingSuccess(token) => Future.successful(Map("Authorization" -> s"Bearer ${token}"))
          case GithubAuthVendingFailure(error) => Future.failed(error)
        }
      }
    }
  }

}
