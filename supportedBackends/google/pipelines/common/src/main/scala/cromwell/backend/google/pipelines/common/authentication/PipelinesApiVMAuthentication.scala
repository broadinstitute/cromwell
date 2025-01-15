package cromwell.backend.google.pipelines.common.authentication

import cats.data.Validated.{Invalid, Valid}
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.DockerCredentials

object PipelinesApiDockerCredentials {

  def apply(dockerCredentials: DockerCredentials, googleConfig: GoogleConfiguration): PipelinesApiDockerCredentials = {
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
        new PipelinesApiDockerCredentials(token = dockerCredentials.token,
                                          keyName = dockerCredentials.keyName,
                                          authName = dockerCredentials.authName
        )
    }
  }
}

/**
 * Authentication information to pull docker images as the user.
 */
case class PipelinesApiDockerCredentials(override val token: String,
                                         override val keyName: Option[String],
                                         override val authName: Option[String]
) extends DockerCredentials(token = token, keyName = keyName, authName = authName)
