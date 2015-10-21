package cromwell.engine.backend.jes

import java.net.URL

import com.typesafe.config.ConfigFactory
import cromwell.engine.backend.jes.authentication.{JesAuthMode, JesDockerCredentials}
import cromwell.util.ConfigUtil._
import cromwell.util.ProductionDockerConfiguration
import cromwell.util.google.{ClientSecrets, SimpleClientSecrets}

import scala.language.postfixOps
import scalaz.Scalaz._
import scalaz.Validation.FlatMap._
import scalaz._

case class JesAttributes(applicationName: String,
                         project: String,
                         executionBucket: String,
                         endpointUrl: URL,
                         authMode: JesAuthMode,
                         dockerCredentials: Option[JesDockerCredentials],
                         googleSecrets: Option[ClientSecrets]) {

  val localizeWithRefreshToken = googleSecrets.isDefined
  val isDockerAuthenticated = dockerCredentials.isDefined
}

object JesAttributes {

  private val jesKeys = Set(
    "applicationName",
    "project",
    "baseExecutionBucket",
    "endpointUrl",
    "maximumPollingInterval",
    "dockerAccount",
    "dockerToken",
    "localizeWithRefreshToken"
  )

  private val googleKeys = Set(
    "authScheme",
    "userAuth",
    "localizeWithRefreshToken.client_id",
    "localizeWithRefreshToken.client_secret",
    "serviceAuth.serviceAccountId",
    "serviceAuth.p12File"
  )

  private val context = "Jes"

  def apply(): JesAttributes = {
    val jesConf = ConfigFactory.load.getConfig("backend").getConfig("jes")
    val googleConf = ConfigFactory.load.getConfig("google")
    val clientSecrets = googleConf.getConfigOption("localizeWithRefreshToken") validateAny {
      _ map { config => SimpleClientSecrets(config.getString("client_id"), config.getString("client_secret")) }
    }

    jesConf.warnNotRecognized(jesKeys, context)
    googleConf.warnNotRecognized(googleKeys, "google")

    val applicationName: ValidationNel[String, String] = jesConf.validateString("applicationName")
    val project: ValidationNel[String, String] = jesConf.validateString("project")
    val executionBucket: ValidationNel[String, String] = jesConf.validateString("baseExecutionBucket")
    val endpointUrl: ValidationNel[String, URL] = jesConf.validateURL("endpointUrl")
    val authMode: ValidationNel[String, JesAuthMode] = googleConf.validateString("authScheme") flatMap {
      _ validateAny JesAuthMode.fromString
    }

    val jesDockerCredentials = ProductionDockerConfiguration.dockerConf map JesDockerCredentials.apply

    (applicationName |@| project |@| executionBucket |@| endpointUrl |@| authMode |@| clientSecrets) {
      JesAttributes(_, _, _, _, _, jesDockerCredentials, _)
    } match {
      case Success(r) => r
      case Failure(f) =>
        val errorMessages = f.list.mkString(", ")
        throw new IllegalArgumentException(s"Jes Configuration is not valid: Errors: $errorMessages")
    }
  }

}