package cromwell.backend.google.pipelines.common

import java.net.URL

import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
import com.typesafe.config.{Config, ConfigValue}
import common.exception.MessageAggregation
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.{BatchRequestTimeoutConfiguration, GcsTransferConfiguration, MemoryRetryConfiguration, VirtualPrivateCloudConfiguration}
import cromwell.backend.google.pipelines.common.authentication.PipelinesApiAuths
import cromwell.backend.google.pipelines.common.callcaching.{CopyCachedOutputs, PipelinesCacheHitDuplicationStrategy, UseOriginalCachedOutputs}
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.{refineMV, refineV}
import net.ceedubs.ficus.Ficus._
import org.slf4j.{Logger, LoggerFactory}

import scala.collection.JavaConverters._
import scala.concurrent.duration._
import scala.util.Try
import scala.util.matching.Regex

case class PipelinesApiConfigurationAttributes(project: String,
                                               computeServiceAccount: String,
                                               auths: PipelinesApiAuths,
                                               restrictMetadataAccess: Boolean,
                                               executionBucket: String,
                                               endpointUrl: URL,
                                               maxPollingInterval: Int,
                                               qps: Int Refined Positive,
                                               cacheHitDuplicationStrategy: PipelinesCacheHitDuplicationStrategy,
                                               requestWorkers: Int Refined Positive,
                                               logFlushPeriod: Option[FiniteDuration],
                                               gcsTransferConfiguration: GcsTransferConfiguration,
                                               virtualPrivateCloudConfiguration: Option[VirtualPrivateCloudConfiguration],
                                               batchRequestTimeoutConfiguration: BatchRequestTimeoutConfiguration,
                                               memoryRetryConfiguration: Option[MemoryRetryConfiguration])

object PipelinesApiConfigurationAttributes {

  /**
    * @param transferAttempts This is the number of attempts, not retries, hence it is positive.
    */
  case class GcsTransferConfiguration(transferAttempts: Int Refined Positive, parallelCompositeUploadThreshold: String)

  final case class VirtualPrivateCloudConfiguration(name: String, subnetwork: Option[String], auth: GoogleAuthMode)
  final case class BatchRequestTimeoutConfiguration(readTimeoutMillis: Option[Int Refined Positive], connectTimeoutMillis: Option[Int Refined Positive])
  final case class MemoryRetryConfiguration(errorKeys: List[String], multiplier: GreaterEqualRefined)


  lazy val Logger = LoggerFactory.getLogger("PipelinesApiConfiguration")

  val GenomicsApiDefaultQps = 1000
  val DefaultGcsTransferAttempts = refineMV[Positive](3)

  lazy val DefaultMemoryRetryFactor: GreaterEqualRefined = refineMV[GreaterEqualOne](2.0)

  private val papiKeys = Set(
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
    "genomics.parallel-composite-upload-threshold",
    "dockerhub",
    "dockerhub.account",
    "dockerhub.token",
    "filesystems",
    "filesystems.gcs.auth",
    "filesystems.gcs.caching.duplication-strategy",
    "concurrent-job-limit",
    "request-workers",
    "batch-requests.timeouts.read",
    "batch-requests.timeouts.connect",
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
    "default-runtime-attributes.zones",
    "virtual-private-cloud",
    "virtual-private-cloud.network-label-key",
    "virtual-private-cloud.subnetwork-label-key",
    "virtual-private-cloud.auth",
    "memory-retry.error-keys",
    "memory-retry.multiplier",
  )

  private val deprecatedJesKeys: Map[String, String] = Map(
    "genomics.default-zones" -> "default-runtime-attributes.zones"
  )

  def apply(googleConfig: GoogleConfiguration, backendConfig: Config, backendName: String): PipelinesApiConfigurationAttributes = {

    def vpcErrorMessage(missingKeys: List[String]) = s"Virtual Private Cloud configuration is invalid. Missing keys: `${missingKeys.mkString(",")}`.".invalidNel

    def validateVPCConfig(networkOption: Option[String], subnetworkOption: Option[String], authOption: Option[String]): ErrorOr[Option[VirtualPrivateCloudConfiguration]] = {
      (networkOption,subnetworkOption, authOption) match {
        case (Some(network), _, Some(auth)) => googleConfig.auth(auth) match {
          case Valid(validAuth) => Option(VirtualPrivateCloudConfiguration(network, subnetworkOption, validAuth)).validNel
          case Invalid(error) => s"Auth $auth is not valid for Virtual Private Cloud configuration. Reason: $error" .invalidNel
        }
        case (Some(_), _, None) => vpcErrorMessage(List("auth"))
        case (None, _, Some(_)) => vpcErrorMessage(List("network-label-key"))
        case (None, Some(_), None) => vpcErrorMessage(List("network-label-key", "auth"))
        case (None, None, None) => None.validNel
      }
    }

    def validateMemoryRetryConfig(errorKeys: Option[List[String]], multiplier: Option[GreaterEqualRefined]): ErrorOr[Option[MemoryRetryConfiguration]] = {
      (errorKeys, multiplier) match {
        case (Some(keys), Some(mul)) => Option(MemoryRetryConfiguration(keys, mul)).validNel
        case (Some(keys), None) => Option(MemoryRetryConfiguration(keys, DefaultMemoryRetryFactor)).validNel
        case (None, Some(_)) => "memory-retry configuration is invalid. No error-keys provided.".invalidNel
        case (None, None) => None.validNel
      }
    }

    val configKeys = backendConfig.entrySet().asScala.toSet map { entry: java.util.Map.Entry[String, ConfigValue] => entry.getKey }
    warnNotRecognized(configKeys, papiKeys, backendName, Logger)

    def warnDeprecated(keys: Set[String], deprecated: Map[String, String], context: String, logger: Logger) = {
      val deprecatedKeys = keys.intersect(deprecated.keySet)
      deprecatedKeys foreach { key => logger.warn(s"Found deprecated configuration key $key, replaced with ${deprecated.get(key)}") }
    }

    warnDeprecated(configKeys, deprecatedJesKeys, backendName, Logger)

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

    val parallelCompositeUploadThreshold = validateGsutilMemorySpecification(backendConfig, "genomics.parallel-composite-upload-threshold")

    val localizationAttempts: ErrorOr[Int Refined Positive] = backendConfig.as[Option[Int]]("genomics.localization-attempts")
      .map(attempts => validatePositiveInt(attempts, "genomics.localization-attempts"))
      .getOrElse(DefaultGcsTransferAttempts.validNel)

    val gcsTransferConfiguration: ErrorOr[GcsTransferConfiguration] =
      (localizationAttempts, parallelCompositeUploadThreshold) mapN GcsTransferConfiguration.apply

    val vpcNetworkLabel: ErrorOr[Option[String]] = validate { backendConfig.getAs[String]("virtual-private-cloud.network-label-key") }
    val vpcSubnetworkLabel: ErrorOr[Option[String]] = validate { backendConfig.getAs[String]("virtual-private-cloud.subnetwork-label-key") }
    val vpcAuth: ErrorOr[Option[String]] = validate { backendConfig.getAs[String]("virtual-private-cloud.auth")}

    val virtualPrivateCloudConfiguration: ErrorOr[Option[VirtualPrivateCloudConfiguration]] = {
      (vpcNetworkLabel, vpcSubnetworkLabel, vpcAuth) flatMapN  validateVPCConfig
    }

    val batchRequestsReadTimeout = readOptionalPositiveMillisecondsIntFromDuration(backendConfig, "batch-requests.timeouts.read")
    val batchRequestsConnectTimeout = readOptionalPositiveMillisecondsIntFromDuration(backendConfig, "batch-requests.timeouts.connect")

    val batchRequestTimeoutConfigurationValidation = (batchRequestsReadTimeout, batchRequestsConnectTimeout) mapN { (read, connect) =>
      BatchRequestTimeoutConfiguration(readTimeoutMillis = read, connectTimeoutMillis = connect)
    }

    val memoryRetryMultiplier: ErrorOr[Option[GreaterEqualRefined]] = validatePositiveOptionDouble(
      backendConfig.as[Option[Double]]("memory-retry.multiplier"),
      configPath = "memory-retry.multiplier"
    )

    val memoryRetryKeys: ErrorOr[Option[List[String]]] = validate { backendConfig.as[Option[List[String]]]("memory-retry.error-keys") }

    val memoryRetryConfig: ErrorOr[Option[MemoryRetryConfiguration]] = {
      (memoryRetryKeys, memoryRetryMultiplier) flatMapN validateMemoryRetryConfig
    }

    def authGoogleConfigForPapiConfigurationAttributes(project: String,
                                                       bucket: String,
                                                       endpointUrl: URL,
                                                       genomicsName: String,
                                                       restrictMetadata: Boolean,
                                                       gcsName: String,
                                                       qps: Int Refined Positive,
                                                       cacheHitDuplicationStrategy: PipelinesCacheHitDuplicationStrategy,
                                                       requestWorkers: Int Refined Positive,
                                                       gcsTransferConfiguration: GcsTransferConfiguration,
                                                       virtualPrivateCloudConfiguration: Option[VirtualPrivateCloudConfiguration],
                                                       batchRequestTimeoutConfiguration: BatchRequestTimeoutConfiguration,
                                                       memoryRetryConfig: Option[MemoryRetryConfiguration]): ErrorOr[PipelinesApiConfigurationAttributes] =
      (googleConfig.auth(genomicsName), googleConfig.auth(gcsName)) mapN {
        (genomicsAuth, gcsAuth) =>
          PipelinesApiConfigurationAttributes(
            project = project,
            computeServiceAccount = computeServiceAccount,
            auths = PipelinesApiAuths(genomicsAuth, gcsAuth),
            restrictMetadataAccess = restrictMetadata,
            executionBucket = bucket,
            endpointUrl = endpointUrl,
            maxPollingInterval = maxPollingInterval,
            qps = qps,
            cacheHitDuplicationStrategy = cacheHitDuplicationStrategy,
            requestWorkers = requestWorkers,
            logFlushPeriod = logFlushPeriod,
            gcsTransferConfiguration = gcsTransferConfiguration,
            virtualPrivateCloudConfiguration = virtualPrivateCloudConfiguration,
            batchRequestTimeoutConfiguration = batchRequestTimeoutConfiguration,
            memoryRetryConfiguration = memoryRetryConfig,
          )
    }

    (project,
      executionBucket,
      endpointUrl,
      genomicsAuthName,
      genomicsRestrictMetadataAccess,
      gcsFilesystemAuthName,
      qpsValidation,
      duplicationStrategy,
      requestWorkers,
      gcsTransferConfiguration,
      virtualPrivateCloudConfiguration,
      batchRequestTimeoutConfigurationValidation,
      memoryRetryConfig,
    ) flatMapN authGoogleConfigForPapiConfigurationAttributes match {
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

  def validateGsutilMemorySpecification(config: Config, configPath: String): ErrorOr[String] = {
    val entry = config.as[Option[String]](configPath)
    entry match {
      case None => "0".validNel
      case Some(v @ GsutilHumanBytes(_, _)) => v.validNel
      case Some(bad) => s"Invalid gsutil memory specification in Cromwell configuration at path '$configPath': '$bad'".invalidNel
    }
  }

  def validatePositiveInt(n: Int, configPath: String) = {
    refineV[Positive](n) match {
      case Left(_) => s"Value $n for $configPath is not strictly positive".invalidNel
      case Right(refined) => refined.validNel
    }
  }

  def readOptionalPositiveMillisecondsIntFromDuration(backendConfig: Config, configPath: String): ErrorOr[Option[Int Refined Positive]] = {

    def validate(n: FiniteDuration) = {
      val result: ErrorOr[Int Refined Positive] = Try(n.toMillis.toInt).toErrorOr flatMap { millisInt =>
        refineV[Positive](millisInt) match {
          case Left(_) => s"Value $n for $configPath is not strictly positive".invalidNel
          case Right(refined) => refined.validNel
        }
      }

      result.contextualizeErrors(s"Parse '$configPath' value $n as a positive Int (in milliseconds)")
    }

    backendConfig.as[Option[FiniteDuration]](configPath) match {
      case Some(value) => validate(value).map(Option.apply)
      case None => None.validNel
    }
  }

  def validatePositiveOptionDouble(value: Option[Double], configPath: String): ErrorOr[Option[GreaterEqualRefined]] = {
    value match {
      case Some(n) => refineV[GreaterEqualOne](n) match {
        case Left(_) => s"Value $n for $configPath should be greater than 1.0.".invalidNel
        case Right(refined) => Option(refined).validNel
      }
      case None => None.validNel
    }
  }

  // Copy/port of gsutil's "_GenerateSuffixRegex"
  private [common] lazy val GsutilHumanBytes: Regex = {
    val _EXP_STRINGS = List(
      List("B", "bit"),
      List("KiB", "Kibit", "K"),
      List("MiB", "Mibit", "M"),
      List("GiB", "Gibit", "G"),
      List("TiB", "Tibit", "T"),
      List("PiB", "Pibit", "P"),
      List("EiB", "Eibit", "E"),
    )

    val suffixes = for {
      unit <- _EXP_STRINGS
      name <- unit
    } yield name

    // Differs from the Python original in a couple of ways:
    //
    // * The Python original uses named groups which are not supported in Scala regexes.
    //   (?P<num>\d*\.\d+|\d+)\s*(?P<suffix>%s)?
    //
    // * The Python original lowercases both the units and the human string before running the matcher.
    //   This Scala version turns on the (?i) case insensitive matching regex option instead.
    val orSuffixes = suffixes.mkString("|")
    "(?i)(\\d*\\.\\d+|\\d+)\\s*(%s)?".format(orSuffixes).r
  }
}
