package cromwell.backend.google.pipelines.common.authentication

import cats.data.Validated.{Invalid, Valid}
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.ClientSecrets
import cromwell.core.DockerCredentials
import spray.json.{JsString, JsValue}
/**
 * Interface for Authentication information that can be included as a json object in the file uploaded to GCS
 * upon workflow creation and used in the VM.
 */
sealed trait PipelinesApiAuthObject {
  def context: String
  def map: Map[String, JsValue]

  def toMap: Map[String, Map[String, JsValue]] =  Map(context -> map)
}

/**
 * Authentication information for data (de)localization as the user.
 */
case class GcsLocalizing(clientSecrets: ClientSecrets, token: String) extends PipelinesApiAuthObject {
  override val context = "boto"
  override val map = Map(
    "client_id" -> JsString(clientSecrets.clientId),
    "client_secret" -> JsString(clientSecrets.clientSecret),
    "refresh_token" -> JsString(token)
  )
}

object PipelinesApiDockerCredentials {

  def apply(dockerCredentials: DockerCredentials, googleConfig: GoogleConfiguration): PipelinesApiDockerCredentials = {
    // If there's an encryption key defined there must be a valid auth defined to encrypt it.
    val authValidation = dockerCredentials.keyName match {
      case None => ().validNel // fine
      case _ =>
        for {
          authName <- dockerCredentials.authName.toErrorOr("KMS Encryption key defined for private Docker but no auth specified")
          _ <- googleConfig.auth(authName)
        } yield ()
    }

    authValidation match {
      case Invalid(errors) =>
        throw new RuntimeException(errors.toList.mkString(", "))
      case Valid(_) =>
        new PipelinesApiDockerCredentials(
          token = dockerCredentials.token,
          keyName = dockerCredentials.keyName,
          authName = dockerCredentials.authName)
    }
  }
}

/**
 * Authentication information to pull docker images as the user.
 */
case class PipelinesApiDockerCredentials(override val token: String,
                                         override val keyName: Option[String],
                                         override val authName: Option[String])
  extends DockerCredentials(token = token, keyName = keyName, authName = authName) with PipelinesApiAuthObject {

  override val context = "docker"
  override val map = Map(
    "token" -> JsString(token)
  )
}
