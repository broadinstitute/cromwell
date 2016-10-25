package cromwell.backend.impl.jes.authentication

import cromwell.core.DockerCredentials
import cromwell.filesystems.gcs.auth.ClientSecrets
import spray.json.{JsString, JsValue}

/**
 * Interface for Authentication information that can be included in the json file uploaded to GCS
 * upon workflow creation and used in the VM.
 */
sealed trait JesAuthInformation {
  def context: String
  def map: Map[String, JsValue]

  def toMap: Map[String, Map[String, JsValue]] =  Map(context -> map)
}

/**
 * Authentication information for data (de)localization as the user.
 */
case class GcsLocalizing(clientSecrets: ClientSecrets, token: String) extends JesAuthInformation {
  override val context = "boto"
  override val map = Map(
    "client_id" -> JsString(clientSecrets.clientId),
    "client_secret" -> JsString(clientSecrets.clientSecret),
    "refresh_token" -> JsString(token)
  )
}

object JesDockerCredentials {
  def apply(dockerCredentials: DockerCredentials): JesDockerCredentials = apply(dockerCredentials.account, dockerCredentials.token)
  def apply(account: String, token: String): JesDockerCredentials = new JesDockerCredentials(account, token)
}

/**
 * Authentication information to pull docker images as the user.
 */
class JesDockerCredentials(account: String, token: String) extends DockerCredentials(account, token) with JesAuthInformation {
  override val context = "docker"
  override val map = Map(
    "account" -> JsString(account),
    "token" -> JsString(token)
  )
}
