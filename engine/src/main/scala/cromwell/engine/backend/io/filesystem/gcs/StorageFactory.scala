package cromwell.engine.backend.io.filesystem.gcs

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.util.DateTime
import com.google.api.services.storage.Storage
import com.google.api.services.storage.model.Bucket.Owner
import cromwell.core.WorkflowOptions
import cromwell.engine.io.gcs.GoogleConfiguration
import cromwell.util.google.GoogleCredentialFactory

import scala.language.postfixOps
import scala.util.Try

object StorageFactory {

  val RefreshTokenOptionKey = "refresh_token"

  lazy val cromwellAuthenticated: Try[Storage] = {
    for {
      cromwellCredentials <- Try(GoogleCredentialFactory.fromCromwellAuthScheme)
      conf <- GoogleConfiguration.gcloudConf
    } yield apply(conf.appName, cromwellCredentials)
  }

  def userAuthenticated(workflowOptions: WorkflowOptions): Try[Storage] = for {
    conf <- GoogleConfiguration.gcloudConf
    token <- workflowOptions.get(RefreshTokenOptionKey)
    userCredentials <- GoogleCredentialFactory.fromUserAuthScheme(token)
  } yield apply(conf.appName, userCredentials)

  private def apply(appName: String, credential: Credential): Storage = {
    new Storage.Builder(GoogleCredentialFactory.httpTransport, GoogleCredentialFactory.jsonFactory, credential).setApplicationName(appName).build()
  }

  case class GcsBucketInfo(bucketName: String, location: String, timeCreated: DateTime, owner: Owner)
}
