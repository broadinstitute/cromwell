package cromwell.services.auth

import akka.actor.ActorRef
import akka.pattern.{AskSupport, AskTimeoutException}
import akka.util.Timeout
import com.typesafe.scalalogging.StrictLogging
import cromwell.languages.util.ImportResolver.{GithubImportAuthProvider, ImportAuthProvider}
import cromwell.services.auth.GithubAuthVending.{
  GithubAuthRequest,
  GithubAuthTokenResponse,
  GithubAuthVendingFailure,
  GithubAuthVendingResponse,
  NoGithubAuthResponse
}

import scala.concurrent.{ExecutionContext, Future}

trait GithubAuthVendingSupport extends AskSupport with StrictLogging {

  def serviceRegistryActor: ActorRef

  implicit val timeout: Timeout
  implicit val ec: ExecutionContext

  def importAuthProvider(token: String): ImportAuthProvider = new GithubImportAuthProvider {
    override def authHeader(): Future[Map[String, String]] =
      serviceRegistryActor
        .ask(GithubAuthRequest(token))
        .mapTo[GithubAuthVendingResponse]
        .recoverWith {
          case e: AskTimeoutException =>
            Future.failed(new Exception(s"Unable to resolve github auth token within allowed time", e))
          case e: Throwable =>
            // This "should" never happen. If it does, let's make it obvious and trigger our alerting:
            logger.error("Programmer error: Unexpected failure to resolve github auth token", e)
            Future.failed(new Exception("Failed to resolve github auth token", e))
        }
        .flatMap {
          case GithubAuthTokenResponse(token) => Future.successful(Map("Authorization" -> s"Bearer ${token}"))
          case NoGithubAuthResponse => Future.successful(Map.empty)
          case GithubAuthVendingFailure(error) =>
            Future.failed(new Exception("Failed to resolve github auth token", error))
        }
  }
}
