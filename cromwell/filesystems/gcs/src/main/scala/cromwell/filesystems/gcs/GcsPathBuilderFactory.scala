package cromwell.filesystems.gcs

import akka.actor.ActorSystem
import com.google.api.gax.retrying.RetrySettings
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration
import com.google.auth.Credentials
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.cloudsupport.gcp.gcs.GcsStorage
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import org.threeten.bp.Duration

import scala.concurrent.ExecutionContext

final case class GcsPathBuilderFactory private(authMode: GoogleAuthMode,
                                               applicationName: String,
                                               retrySettings: RetrySettings,
                                               cloudStorageConfiguration: CloudStorageConfiguration)
  extends PathBuilderFactory {

  def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext) = {
    GcsPathBuilder.fromAuthMode(authMode, applicationName, retrySettings, GcsStorage.DefaultCloudStorageConfiguration, options)
  }

  /**
    * Ignores the authMode and creates a GcsPathBuilder using the passed credentials directly.
    * Can be used when the Credentials are already available.
    */
  def fromCredentials(options: WorkflowOptions, credentials: Credentials) = {
    GcsPathBuilder.fromCredentials(credentials, applicationName, retrySettings, GcsStorage.DefaultCloudStorageConfiguration, options)
  }
}

object GcsPathBuilderFactory {
  def apply(authMode: GoogleAuthMode,
            applicationName: String,
            retrySettings: Option[RetrySettings]): GcsPathBuilderFactory = {

    val actualRetrySettings: RetrySettings = retrySettings.getOrElse(DefaultRetrySettings)
    val actualCloudStorage = DefaultCloudStorageConfiguration

    new GcsPathBuilderFactory(authMode, applicationName, actualRetrySettings, actualCloudStorage)
  }

  lazy val DefaultRetrySettings: RetrySettings = RetrySettings.newBuilder()
    .setTotalTimeout(Duration.ofSeconds(30))
    .setInitialRetryDelay(Duration.ofMillis(100))
    .setRetryDelayMultiplier(1.1)
    .setMaxRetryDelay(Duration.ofSeconds(1))
    .setJittered(true)
    .setInitialRpcTimeout(Duration.ofMillis(100))
    .setRpcTimeoutMultiplier(1.1)
    .setMaxRpcTimeout(Duration.ofSeconds(5))
    .build()

  lazy val DefaultCloudStorageConfiguration = GcsStorage.DefaultCloudStorageConfiguration
}
