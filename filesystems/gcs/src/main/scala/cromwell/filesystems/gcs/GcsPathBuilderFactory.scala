package cromwell.filesystems.gcs

import akka.actor.ActorSystem
import com.google.api.gax.retrying.RetrySettings
import com.google.common.cache.CacheBuilder
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.cloudsupport.gcp.gcs.GcsStorage
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import cromwell.filesystems.gcs.GcsPathBuilderFactory.DefaultRetrySettings
import cromwell.filesystems.gcs.cache.GcsBucketCache.GcsBucketInformation
import cromwell.filesystems.gcs.cache.GcsBucketInformationPolicies
import cromwell.filesystems.gcs.cache.GcsBucketInformationPolicies.GcsBucketInformationPolicy
import org.threeten.bp.Duration

import scala.concurrent.ExecutionContext
import scala.concurrent.duration.{Duration => SDuration, _}

final case class GcsPathBuilderFactory(globalConfig: Config, instanceConfig: Config)
  extends PathBuilderFactory {
  import net.ceedubs.ficus.Ficus._
  // Parse the configuration and create a GoogleConfiguration
  val googleConf: GoogleConfiguration = GoogleConfiguration(globalConfig)
  // Extract the specified authentication mode for gcs filesystem, if any
  val authModeAsString: String = instanceConfig.as[String]("auth")
  // Validate it against the google configuration
  val authModeValidation: ErrorOr[GoogleAuthMode] = googleConf.auth(authModeAsString)
  val applicationName = googleConf.applicationName

  val authMode = authModeValidation.unsafe(s"Failed to create authentication mode for $authModeAsString")

  val requesterPaysEnabled = instanceConfig.as[Option[Config]]("requester-pays")
  val defaultProject = requesterPaysEnabled.flatMap(_.as[Option[String]]("project"))
  val requesterPaysCacheTTL = requesterPaysEnabled.flatMap(_.as[Option[FiniteDuration]]("cache-ttl")).getOrElse(20.minutes)
  
  def cachePolicy(ttl: FiniteDuration) = {
    val guavaCache = CacheBuilder
      .newBuilder()
      .expireAfterWrite(ttl.length, ttl.unit)
      .build[String, GcsBucketInformation]()

    GcsBucketInformationPolicies.CacheInformation(guavaCache)
  }
  
  // Determines the behavior of bucket information retrieval and caching
  val gcsBucketInformationPolicy: GcsBucketInformationPolicy = requesterPaysCacheTTL match {
    // Default to no lookup at all if there is no config.
    // This preserves the behavior on V1, however it forces to set a config value in V2
    case _ if requesterPaysEnabled.isEmpty => GcsBucketInformationPolicies.Default
    case inf if !inf.isFinite() => GcsBucketInformationPolicies.Default
    case zero if zero == SDuration.Zero => GcsBucketInformationPolicies.NoCache
    case other => cachePolicy(other)
  }

  def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext) = {
    GcsPathBuilder.fromAuthMode(authMode, applicationName, DefaultRetrySettings, GcsStorage.DefaultCloudStorageConfiguration, options, defaultProject, gcsBucketInformationPolicy)
  }
}

object GcsPathBuilderFactory {
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
