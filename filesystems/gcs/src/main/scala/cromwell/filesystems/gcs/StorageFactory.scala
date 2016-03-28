package cromwell.filesystems.gcs

import com.google.api.client.util.DateTime
import com.google.api.services.storage.Storage
import com.google.api.services.storage.model.Bucket.Owner

import scala.language.postfixOps

object StorageFactory {

  def apply(gcsConf: GoogleConfiguration, refreshToken: Option[RefreshToken] = None): Storage = {
    val credential = GoogleCredentialFactory(gcsConf.authMode, GcsScopes, refreshToken)
    new Storage.Builder(credential.getTransport, credential.getJsonFactory, credential).setApplicationName(gcsConf.appName).build()
  }

  case class GcsBucketInfo(bucketName: String, location: String, timeCreated: DateTime, owner: Owner)
}
