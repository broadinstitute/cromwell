package cromwell.services.auth

import akka.actor.ActorRef
import akka.pattern.{AskSupport, AskTimeoutException}
import akka.util.Timeout
import cats.implicits.catsSyntaxValidatedId
import com.typesafe.config.Config
import com.typesafe.scalalogging.StrictLogging
import common.validation.ErrorOr.ErrorOr
import cromwell.languages.util.ImportResolver.{GithubImportAuthProvider, ImportAuthProvider}
import cromwell.services.auth.GithubAuthVending._

import scala.concurrent.{ExecutionContext, Future}

trait GithubAuthVendingSupport extends AskSupport with StrictLogging {

  def serviceRegistryActor: ActorRef

  implicit val ec: ExecutionContext

  def importAuthProvider(token: String)(implicit timeout: Timeout): ImportAuthProvider = new GithubImportAuthProvider {
    override def authHeader(): Future[Map[String, String]] =
      serviceRegistryActor
        .ask(GithubAuthRequest(TerraToken(token)))
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
          case GithubAuthTokenResponse(githubToken) =>
            Future.successful(Map("Authorization" -> s"Bearer ${githubToken.value}"))
          case NoGithubAuthResponse => Future.successful(Map.empty)
          case GithubAuthVendingFailure(error) =>
            Future.failed(new Exception(s"Failed to resolve GitHub auth token. Error: $error"))
        }
  }

  def importAuthProvider(config: Config): ErrorOr[Option[ImportAuthProvider]] =
    // TODO add GCP auth provider
    None.validNel
}
