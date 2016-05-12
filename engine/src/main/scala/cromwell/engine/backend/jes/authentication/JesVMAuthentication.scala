package cromwell.engine.backend.jes.authentication

import cromwell.core.DockerCredentials
import cromwell.filesystems.gcs.ClientSecrets
import spray.json.{JsString, JsValue}

/**
 * Interface for Authentication information that can be included in the json file uploaded to GCS
 * upon workflow creation and used in the VM.
 */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
sealed trait JesAuthInformation {
  def context: String
  def map: Map[String, JsValue]

  def toMap: Map[String, Map[String, JsValue]] =  Map(context -> map)
}

/**
 * Authentication information for data (de)localization as the user.
 */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class GcsLocalizing(clientSecrets: ClientSecrets, token: String) extends JesAuthInformation {
  override val context = "boto"
  override val map = Map(
    "client_id" -> JsString(clientSecrets.clientId),
    "client_secret" -> JsString(clientSecrets.clientSecret),
    "refresh_token" -> JsString(token)
  )
}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object JesDockerCredentials {
  def apply(dockerCredentials: DockerCredentials): JesDockerCredentials = apply(dockerCredentials.account, dockerCredentials.token)
  def apply(account: String, token: String): JesDockerCredentials = new JesDockerCredentials(account, token)
}

/**
 * Authentication information to pull docker images as the user.
 */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
class JesDockerCredentials(account: String, token: String) extends DockerCredentials(account, token) with JesAuthInformation {
  override val context = "docker"
  override val map = Map(
    "account" -> JsString(account),
    "token" -> JsString(token)
  )
}
