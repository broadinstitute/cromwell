package cromwell.backend.google.pipelines.common

import java.net.URL

import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
import com.typesafe.config.{Config, ConfigValue}
import common.exception.MessageAggregation
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.{LocalizationConfiguration, BatchRequestTimeoutConfiguration, VirtualPrivateCloudConfiguration}
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
                                               localizationConfiguration: LocalizationConfiguration,
                                               virtualPrivateCloudConfiguration: Option[VirtualPrivateCloudConfiguration],
                                               batchRequestTimeoutConfiguration: BatchRequestTimeoutConfiguration)

object PipelinesApiConfigurationAttributes {

  /**
    * @param localizationAttempts Also used for de-localization. This is the number of attempts, not retries,
    *                             hence it is positive.
    */
  case class LocalizationConfiguration(localizationAttempts: Int Refined Positive)

  final case class VirtualPrivateCloudConfiguration(name: String, subnetwork: Option[String], auth: GoogleAuthMode)
  final case class BatchRequestTimeoutConfiguration(readTimeoutMillis: Option[Int Refined Positive], connectTimeoutMillis: Option[Int Refined Positive])


  lazy val Logger = LoggerFactory.getLogger("PipelinesApiConfiguration")

  val GenomicsApiDefaultQps = 1000
  val DefaultLocalizationAttempts = refineMV[Positive](3)

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
    "virtual-private-cloud.auth"
  )

  private val deprecatedJesKeys: Map[String, String] = Map(
    "genomics.default-zones" -> "default-runtime-attributes.zones"
  )

  private val context = "Jes"

  def apply(googleConfig: GoogleConfiguration, backendConfig: Config): PipelinesApiConfigurationAttributes = {

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

    val configKeys = backendConfig.entrySet().asScala.toSet map { entry: java.util.Map.Entry[String, ConfigValue] => entry.getKey }
    warnNotRecognized(configKeys, papiKeys, context, Logger)

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

    def authGoogleConfigForPapiConfigurationAttributes(project: String,
                                                       bucket: String,
                                                       endpointUrl: URL,
                                                       genomicsName: String,
                                                       restrictMetadata: Boolean,
                                                       gcsName: String,
                                                       qps: Int Refined Positive,
                                                       cacheHitDuplicationStrategy: PipelinesCacheHitDuplicationStrategy,
                                                       requestWorkers: Int Refined Positive,
                                                       localizationConfiguration: LocalizationConfiguration,
                                                       virtualPrivateCloudConfiguration: Option[VirtualPrivateCloudConfiguration],
                                                       batchRequestTimeoutConfiguration: BatchRequestTimeoutConfiguration): ErrorOr[PipelinesApiConfigurationAttributes] =
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
            localizationConfiguration = localizationConfiguration,
            virtualPrivateCloudConfiguration = virtualPrivateCloudConfiguration,
            batchRequestTimeoutConfiguration = batchRequestTimeoutConfiguration
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
      localizationConfiguration,
      virtualPrivateCloudConfiguration,
      batchRequestTimeoutConfigurationValidation
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
}
