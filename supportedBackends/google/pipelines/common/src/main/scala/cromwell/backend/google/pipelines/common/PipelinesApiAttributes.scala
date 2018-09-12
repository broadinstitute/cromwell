package cromwell.backend.google.pipelines.common

import java.net.URL

import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
import com.typesafe.config.{Config, ConfigValue}
import common.exception.MessageAggregation
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.backend.google.pipelines.common.PipelinesApiAttributes.LocalizationConfiguration
import cromwell.backend.google.pipelines.common.authentication.PipelinesApiAuths
import cromwell.backend.google.pipelines.common.callcaching.{CopyCachedOutputs, PipelinesCacheHitDuplicationStrategy, UseOriginalCachedOutputs}
import cromwell.cloudsupport.gcp.GoogleConfiguration
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.{refineMV, refineV}
import net.ceedubs.ficus.Ficus._
import org.slf4j.{Logger, LoggerFactory}

import scala.collection.JavaConverters._
import scala.concurrent.duration._

case class PipelinesApiAttributes(project: String,
                                  computeServiceAccount: String,
                                  auths: PipelinesApiAuths,
                                  restrictMetadataAccess: Boolean,
                                  executionBucket: String,
                                  endpointUrl: URL,
                                  maxPollingInterval: Int,
                                  qps: Int Refined Positive,
                                  duplicationStrategy: PipelinesCacheHitDuplicationStrategy,
                                  requestWorkers: Int Refined Positive,
                                  logFlushPeriod: Option[FiniteDuration],
                                  localizationConfiguration: LocalizationConfiguration)

object PipelinesApiAttributes {

  /**
    * @param localizationAttempts Also used for de-localization. This is the number of attempts, not retries,
    *                             hence it is positive.
    */
  case class LocalizationConfiguration(localizationAttempts: Int Refined Positive)

  lazy val Logger = LoggerFactory.getLogger("JesAttributes")

  val GenomicsApiDefaultQps = 1000
  val DefaultLocalizationAttempts = refineMV[Positive](3)

  private val jesKeys = Set(
    "project",
    "root",
    "maximum-polling-interval",
    "genomics",
    "genomics.compute-service-account",
    "genomics.auth",
    "genomics.restrict-metadata-access",
    "genomics.endpoint-url",
    "genomics-api-queries-per-100-seconds",
    "genomics.localization-attempts",
    "dockerhub",
    "dockerhub.account",
    "dockerhub.token",
    "filesystems",
    "filesystems.gcs.auth",
    "filesystems.gcs.caching.duplication-strategy",
    "concurrent-job-limit",
    "request-workers",
    "default-runtime-attributes",
    "default-runtime-attributes.cpu",
    "default-runtime-attributes.failOnStderr",
    "default-runtime-attributes.continueOnReturnCode",
    "default-runtime-attributes.docker",
    "default-runtime-attributes.memory",
    "default-runtime-attributes.bootDiskSizeGb",
    "default-runtime-attributes.disks",
    "default-runtime-attributes.noAddress",
    "default-runtime-attributes.preemptible",
    "default-runtime-attributes.zones"
  )

  private val deprecatedJesKeys: Map[String, String] = Map(
    "genomics.default-zones" -> "default-runtime-attributes.zones"
  )

  private val context = "Jes"

  def apply(googleConfig: GoogleConfiguration, backendConfig: Config): PipelinesApiAttributes = {
    val configKeys = backendConfig.entrySet().asScala.toSet map { entry: java.util.Map.Entry[String, ConfigValue] => entry.getKey }
    warnNotRecognized(configKeys, jesKeys, context, Logger)

    def warnDeprecated(keys: Set[String], deprecated: Map[String, String], context: String, logger: Logger) = {
      val deprecatedKeys = keys.intersect(deprecated.keySet)
      deprecatedKeys foreach { key => logger.warn(s"Found deprecated configuration key $key, replaced with ${deprecated.get(key)}") }
    }

    warnDeprecated(configKeys, deprecatedJesKeys, context, Logger)

    val project: ErrorOr[String] = validate { backendConfig.as[String]("project") }
    val executionBucket: ErrorOr[String] = validate { backendConfig.as[String]("root") }
    val endpointUrl: ErrorOr[URL] = validate { backendConfig.as[URL]("genomics.endpoint-url") }
    val maxPollingInterval: Int = backendConfig.as[Option[Int]]("maximum-polling-interval").getOrElse(600)
    val computeServiceAccount: String = backendConfig.as[Option[String]]("genomics.compute-service-account").getOrElse("default")
    val genomicsAuthName: ErrorOr[String] = validate { backendConfig.as[String]("genomics.auth") }
    val genomicsRestrictMetadataAccess: ErrorOr[Boolean] = validate { backendConfig.as[Option[Boolean]]("genomics.restrict-metadata-access").getOrElse(false) }
    val gcsFilesystemAuthName: ErrorOr[String] = validate { backendConfig.as[String]("filesystems.gcs.auth") }
    val qpsValidation = validateQps(backendConfig)
    val duplicationStrategy = validate { backendConfig.as[Option[String]]("filesystems.gcs.caching.duplication-strategy").getOrElse("copy") match {
      case "copy" => CopyCachedOutputs
      case "reference" => UseOriginalCachedOutputs
      case other => throw new IllegalArgumentException(s"Unrecognized caching duplication strategy: $other. Supported strategies are copy and reference. See reference.conf for more details.")
    } }
    val requestWorkers: ErrorOr[Int Refined Positive] = validatePositiveInt(backendConfig.as[Option[Int]]("request-workers").getOrElse(3), "request-workers")
    val logFlushPeriod: Option[FiniteDuration] = backendConfig.as[Option[FiniteDuration]]("log-flush-period") match {
      case Some(duration) if duration.isFinite() => Option(duration)
      // "Inf" disables upload
      case Some(_) => None
      // Defaults to 1 minute
      case None => Option(1.minute)
    }

    val localizationConfiguration: ErrorOr[LocalizationConfiguration] =
      backendConfig.as[Option[Int]]("genomics.localization-attempts")
        .map(attempts => validatePositiveInt(attempts, "genomics.localization-attempts"))
        .map(_.map(LocalizationConfiguration.apply))
        .getOrElse(LocalizationConfiguration(DefaultLocalizationAttempts).validNel)


    def authGoogleConfigForJesAttributes(project: String,
                                         bucket: String,
                                         endpointUrl: URL,
                                         genomicsName: String,
                                         restrictMetadata: Boolean,
                                         gcsName: String,
                                         qps: Int Refined Positive,
                                         cachingStrategy: PipelinesCacheHitDuplicationStrategy,
                                         requestWorkers: Int Refined Positive,
                                         localizationConfiguration: LocalizationConfiguration): ErrorOr[PipelinesApiAttributes] = (googleConfig.auth(genomicsName), googleConfig.auth(gcsName)) mapN {
      (genomicsAuth, gcsAuth) => PipelinesApiAttributes(project, computeServiceAccount, PipelinesApiAuths(genomicsAuth, gcsAuth), restrictMetadata, bucket, endpointUrl, maxPollingInterval, qps, cachingStrategy, requestWorkers, logFlushPeriod, localizationConfiguration)
    }

    (project, executionBucket, endpointUrl, genomicsAuthName, genomicsRestrictMetadataAccess, gcsFilesystemAuthName,
      qpsValidation, duplicationStrategy, requestWorkers, localizationConfiguration) flatMapN authGoogleConfigForJesAttributes match {
      case Valid(r) => r
      case Invalid(f) =>
        throw new IllegalArgumentException with MessageAggregation {
          override val exceptionContext = "Google Pipelines API configuration is not valid: Errors"
          override val errorMessages = f.toList
        }
    }
  }

  def validateQps(config: Config): ErrorOr[Int Refined Positive] = {
    import eu.timepit.refined._

    val qp100s = config.as[Option[Int]]("genomics-api-queries-per-100-seconds").getOrElse(GenomicsApiDefaultQps)
    val qpsCandidate = qp100s / 100

    refineV[Positive](qpsCandidate) match {
      case Left(_) => s"Calculated QPS for Google Genomics API ($qpsCandidate/s) was not a positive integer (supplied value was $qp100s per 100s)".invalidNel
      case Right(refined) => refined.validNel
    }
  }

  def validatePositiveInt(n: Int, configPath: String) = {
    refineV[Positive](n) match {
      case Left(_) => s"Value $n for $configPath is not strictly positive".invalidNel
      case Right(refined) => refined.validNel
    }
  }
}
