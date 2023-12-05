package cromwell.backend.google.batch.authentication

import cats.data.Validated.{Invalid, Valid}
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.DockerCredentials
import spray.json.{JsString, JsValue}

/**
  * Interface for Authentication information that can be included as a json object in the file uploaded to GCS
  * upon workflow creation and used in the VM.
  */
sealed trait GcpBatchAuthObject {
  def context: String

  def map: Map[String, JsValue]

  def toMap: Map[String, Map[String, JsValue]] = Map(context -> map)
}

object GcpBatchDockerCredentials {

  def apply(dockerCredentials: DockerCredentials, googleConfig: GoogleConfiguration): GcpBatchDockerCredentials = {
    // If there's an encryption key defined there must be a valid auth defined to encrypt it.
    val authValidation = dockerCredentials.keyName match {
      case None => ().validNel // fine
      case _ =>
        for {
          authName <- dockerCredentials.authName.toErrorOr(
            "KMS Encryption key defined for private Docker but no auth specified"
          )
          _ <- googleConfig.auth(authName)
        } yield ()
    }

    authValidation match {
      case Invalid(errors) =>
        throw new RuntimeException(errors.toList.mkString(", "))
      case Valid(_) =>
        new GcpBatchDockerCredentials(token = dockerCredentials.token,
                                      keyName = dockerCredentials.keyName,
                                      authName = dockerCredentials.authName
        )
    }
  }
}

/**
  * Authentication information to pull docker images as the user.
  */
case class GcpBatchDockerCredentials(override val token: String,
                                     override val keyName: Option[String],
                                     override val authName: Option[String]
) extends DockerCredentials(token = token, keyName = keyName, authName = authName)
    with GcpBatchAuthObject {

  override val context = "docker"
  override val map = Map(
    "token" -> JsString(token)
  )
}
