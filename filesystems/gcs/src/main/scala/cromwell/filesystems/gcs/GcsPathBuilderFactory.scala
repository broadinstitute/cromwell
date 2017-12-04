package cromwell.filesystems.gcs

import akka.actor.ActorSystem
import com.google.api.gax.retrying.RetrySettings
import com.google.auth.Credentials
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.cloudsupport.gcp.gcs.GcsStorage
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory

import scala.concurrent.ExecutionContext

final case class GcsPathBuilderFactory(authMode: GoogleAuthMode,
                                       applicationName: String,
                                       retrySettings: Option[RetrySettings] = None,
                                       cloudStorageConfiguration: CloudStorageConfiguration = GcsStorage.DefaultCloudStorageConfiguration)
  extends PathBuilderFactory {

  def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext) = {
    GcsPathBuilder.fromAuthMode(authMode, applicationName, retrySettings, cloudStorageConfiguration, options)
  }

  /**
    * Ignores the authMode and creates a GcsPathBuilder using the passed credentials directly.
    * Can be used when the Credentials are already available.
    */
  def fromCredentials(options: WorkflowOptions, credentials: Credentials) = {
    GcsPathBuilder.fromCredentials(credentials, applicationName, retrySettings, cloudStorageConfiguration, options)
  }
}
