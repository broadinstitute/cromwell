package cromwell.filesystems.gcs

import akka.actor.ActorSystem
import com.google.api.gax.retrying.RetrySettings
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.cloudsupport.gcp.gcs.GcsStorage
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import org.threeten.bp.Duration

import scala.concurrent.{ExecutionContext, Future}

final case class GcsPathBuilderFactory(globalConfig: Config, instanceConfig: Config)
  extends PathBuilderFactory {
  import net.ceedubs.ficus.Ficus._
  // Parse the configuration and create a GoogleConfiguration
  val googleConf: GoogleConfiguration = GoogleConfiguration(globalConfig)
  // Extract the specified authentication mode for engine gcs filesystem, if any
  val authModeAsString: String = instanceConfig.as[String]("auth")
  // Validate it against the google configuration
  val authModeValidation: ErrorOr[GoogleAuthMode] = googleConf.auth(authModeAsString)
  val applicationName = googleConf.applicationName
  val maxAttempts = instanceConfig.getOrElse("max-attempts", 0)

  val authMode = authModeValidation.unsafe(s"Failed to create authentication mode for $authModeAsString")

  val defaultProject = instanceConfig.as[Option[String]]("project")

  lazy val defaultRetrySettings: RetrySettings = {
    RetrySettings.newBuilder()
      .setMaxAttempts(maxAttempts)
      .setTotalTimeout(Duration.ofSeconds(30))
      .setInitialRetryDelay(Duration.ofMillis(100))
      .setRetryDelayMultiplier(1.1)
      .setMaxRetryDelay(Duration.ofSeconds(1))
      .setInitialRpcTimeout(Duration.ofMillis(100))
      .setRpcTimeoutMultiplier(1.1)
      .setMaxRpcTimeout(Duration.ofSeconds(5))
      .build()
  }

  def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[GcsPathBuilder] = {
    GcsPathBuilder.fromAuthMode(
      authMode,
      applicationName,
      defaultRetrySettings,
      GcsStorage.DefaultCloudStorageConfiguration,
      options,
      defaultProject
    )
  }
}

object GcsPathBuilderFactory {
  lazy val DefaultCloudStorageConfiguration = GcsStorage.DefaultCloudStorageConfiguration
}
