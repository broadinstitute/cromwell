package cromwell.services.auth

import akka.actor.ActorRef
import akka.pattern.AskSupport
import akka.util.Timeout
import cromwell.languages.util.ImportResolver.{GithubImportAuthProvider, ImportAuthProvider}
import cromwell.services.auth.GithubAuthVending.{GithubAuthRequest, GithubAuthTokenResponse, GithubAuthVendingFailure, GithubAuthVendingResponse, NoGithubAuthResponse}

import scala.concurrent.{ExecutionContext, Future}

trait GithubAuthVendingSupport extends AskSupport {

  def serviceRegistryActor: ActorRef

  implicit val timeout: Timeout
  implicit val ec: ExecutionContext

  def importAuthProvider(token: String): ImportAuthProvider = new GithubImportAuthProvider {
    override def authHeader(): Future[Map[String, String]] = {
      serviceRegistryActor.ask(GithubAuthRequest(token)).mapTo[GithubAuthVendingResponse].flatMap {
        case GithubAuthTokenResponse(token) => Future.successful(Map("Authorization" -> s"Bearer ${token}"))
        case NoGithubAuthResponse => Future.successful(Map.empty)
        case GithubAuthVendingFailure(error) => Future.failed(error)
      }
    }
  }
}
