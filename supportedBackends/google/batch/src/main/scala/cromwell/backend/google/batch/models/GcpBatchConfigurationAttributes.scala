package cromwell.backend.google.batch.models

import cats.data.Validated._
import cats.data.{NonEmptyList, Validated}
import cats.implicits._
import com.typesafe.config.{Config, ConfigValue}
import com.typesafe.scalalogging.StrictLogging
import common.exception.MessageAggregation
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.backend.CommonBackendConfigurationAttributes
import cromwell.backend.google.batch._
import cromwell.backend.google.batch.authentication.GcpBatchAuths
import cromwell.backend.google.batch.callcaching.{
  BatchCacheHitDuplicationStrategy,
  CopyCachedOutputs,
  UseOriginalCachedOutputs
}
import cromwell.backend.google.batch.io.GcpBatchReferenceFilesDisk
import cromwell.backend.google.batch.models.GcpBatchConfigurationAttributes.{
  BatchRequestTimeoutConfiguration,
  GcsTransferConfiguration,
  VirtualPrivateCloudConfiguration
}
import cromwell.backend.google.batch.util.GcpBatchReferenceFilesMappingOperations
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.docker.DockerMirroring
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.{refineMV, refineV}
import net.ceedubs.ficus.Ficus._
import org.slf4j.{Logger, LoggerFactory}

import scala.concurrent.duration._
import scala.jdk.CollectionConverters._
import scala.util.matching.Regex
import scala.util.{Failure, Success, Try}

case class GcpBatchConfigurationAttributes(
  project: String,
  computeServiceAccount: String,
  auths: GcpBatchAuths,
  dockerhubToken: String,
  enableFuse: Boolean,
  executionBucket: String,
  location: String,
  maxPollingInterval: Int,
  qps: Int Refined Positive,
  cacheHitDuplicationStrategy: BatchCacheHitDuplicationStrategy,
  requestWorkers: Int Refined Positive,
  batchTimeout: FiniteDuration,
  dockerMirroringOpt: Option[DockerMirroring],
  logFlushPeriod: Option[FiniteDuration],
  gcsTransferConfiguration: GcsTransferConfiguration,
  virtualPrivateCloudConfiguration: VirtualPrivateCloudConfiguration,
  batchRequestTimeoutConfiguration: BatchRequestTimeoutConfiguration,
  referenceFileToDiskImageMappingOpt: Option[Map[String, GcpBatchReferenceFilesDisk]],
  checkpointingInterval: FiniteDuration,
  logsPolicy: GcpBatchLogsPolicy,
  maxTransientErrorRetries: Int
)

object GcpBatchConfigurationAttributes extends GcpBatchReferenceFilesMappingOperations with StrictLogging {

  /**
    * param transferAttempts This is the number of attempts, not retries, hence it is positive.
    */
  case class GcsTransferConfiguration(transferAttempts: Int Refined Positive, parallelCompositeUploadThreshold: String)

  final case class VirtualPrivateCloudLabels(network: String, subnetwork: Option[String], auth: GoogleAuthMode)

  final case class VirtualPrivateCloudLiterals(network: String, subnetwork: Option[String])

  final case class VirtualPrivateCloudConfiguration(labelsOption: Option[VirtualPrivateCloudLabels],
                                                    literalsOption: Option[VirtualPrivateCloudLiterals]
  )

  final case class BatchRequestTimeoutConfiguration(readTimeoutMillis: Option[Int Refined Positive],
                                                    connectTimeoutMillis: Option[Int Refined Positive]
  )

  lazy val Logger: Logger = LoggerFactory.getLogger("BatchConfiguration")

  // Default queries per 100 seconds
  // Batch allows up to 6000 queries per minute (100 requests per second).
  // See: https://cloud.google.com/batch/quotas
  val BatchApiDefaultQp100s = 10000
  val DefaultGcsTransferAttempts: Refined[Int, Positive] = refineMV[Positive](3)

  val checkpointingIntervalKey = "checkpointing-interval"

  private val batchKeys = CommonBackendConfigurationAttributes.commonValidConfigurationAttributeKeys ++ Set(
    "project",
    "root",
    "maximum-polling-interval",
    "batch",
    "batch.location",
    "batch.compute-service-account",
    "batch.auth",
    "batch.restrict-metadata-access",
    "batch.enable-fuse",
    "batch-queries-per-100-seconds",
    "batch.localization-attempts",
    "batch.parallel-composite-upload-threshold",
    "batch.logs-policy",
    "filesystems",
    "filesystems.drs.auth",
    "filesystems.gcs.auth",
    "filesystems.gcs.project",
    "filesystems.gcs.caching.duplication-strategy",
    "concurrent-job-limit",
    "request-workers",
    "batch-timeout",
    "max-transient-error-retries",
    "docker-mirror.dockerhub.enabled",
    "docker-mirror.dockerhub.address",
    "batch-requests.timeouts.read",
    "batch-requests.timeouts.connect",
    "default-runtime-attributes.bootDiskSizeGb",
    "default-runtime-attributes.noAddress",
    "default-runtime-attributes.preemptible",
    "default-runtime-attributes.zones",
    "virtual-private-cloud",
    "virtual-private-cloud.network-name",
    "virtual-private-cloud.subnetwork-name",
    "virtual-private-cloud.network-label-key",
    "virtual-private-cloud.subnetwork-label-key",
    "virtual-private-cloud.auth",
    "reference-disk-localization-manifests",
    "reference-disk-localization.validation.auth",
    checkpointingIntervalKey
  )

  private val deprecatedBatchKeys: Map[String, String] = Map(
    "batch.default-zones" -> "default-runtime-attributes.zones"
  )

  def apply(googleConfig: GoogleConfiguration,
            backendConfig: Config,
            backendName: String
  ): GcpBatchConfigurationAttributes = {

    def vpcErrorMessage(missingKeys: List[String]) =
      s"Virtual Private Cloud configuration is invalid. Missing keys: `${missingKeys.mkString(",")}`.".invalidNel

    def validateVPCLabelsConfig(networkOption: Option[String],
                                subnetworkOption: Option[String],
                                authOption: Option[String]
    ): ErrorOr[Option[VirtualPrivateCloudLabels]] =
      (networkOption, subnetworkOption, authOption) match {
        case (Some(network), _, Some(auth)) =>
          googleConfig.auth(auth) match {
            case Valid(validAuth) =>
              Option(VirtualPrivateCloudLabels(network, subnetworkOption, validAuth)).validNel
            case Invalid(error) =>
              s"Auth $auth is not valid for Virtual Private Cloud configuration. Reason: $error".invalidNel
          }
        case (Some(_), _, None) => vpcErrorMessage(List("auth"))
        case (None, _, Some(_)) => vpcErrorMessage(List("network-label-key"))
        case (None, Some(_), None) => vpcErrorMessage(List("network-label-key", "auth"))
        case (None, None, None) => None.validNel
      }

    def validateVPCLiteralsConfig(networkNameOption: Option[String],
                                  subnetworkNameOption: Option[String]
    ): ErrorOr[Option[VirtualPrivateCloudLiterals]] =
      (networkNameOption, subnetworkNameOption) match {
        case (None, Some(_)) => vpcErrorMessage(List("network-name"))
        case (Some(networkName), _) => Option(VirtualPrivateCloudLiterals(networkName, subnetworkNameOption)).valid
        case (None, None) => None.valid
      }

    def validateVPCConfig(networkNameOption: Option[String],
                          subnetworkNameOption: Option[String],
                          networkLabelOption: Option[String],
                          subnetworkLabelOption: Option[String],
                          authOption: Option[String]
    ): ErrorOr[VirtualPrivateCloudConfiguration] = {
      val vpcLabelsValidation =
        validateVPCLabelsConfig(networkLabelOption, subnetworkLabelOption, authOption)
      val vpcLiteralsValidation =
        validateVPCLiteralsConfig(networkNameOption, subnetworkNameOption)
      (vpcLabelsValidation, vpcLiteralsValidation) mapN VirtualPrivateCloudConfiguration
    }

    val configKeys = backendConfig.entrySet().asScala.toSet map { entry: java.util.Map.Entry[String, ConfigValue] =>
      entry.getKey
    }
    warnNotRecognized(configKeys, batchKeys, backendName, Logger)

    def warnDeprecated(keys: Set[String], deprecated: Map[String, String], context: String, logger: Logger): Unit = {
      val deprecatedKeys = keys.intersect(deprecated.keySet)
      deprecatedKeys foreach { key =>
        logger.warn(s"Found deprecated configuration key $key in $context, replaced with ${deprecated.get(key)}")
      }
    }

    warnDeprecated(configKeys, deprecatedBatchKeys, backendName, Logger)

    val project: ErrorOr[String] = validate {
      backendConfig.as[String]("project")
    }
    val executionBucket: ErrorOr[String] = validate {
      backendConfig.as[String]("root")
    }
    val location: ErrorOr[String] = validate {
      backendConfig.as[String]("batch.location")
    }
    val maxPollingInterval: Int = backendConfig.as[Option[Int]]("maximum-polling-interval").getOrElse(600)
    val computeServiceAccount: String =
      backendConfig.as[Option[String]]("batch.compute-service-account").getOrElse("default")
    val batchAuthName: ErrorOr[String] = validate {
      backendConfig.as[String]("batch.auth")
    }
    val batchRestrictMetadataAccess: ErrorOr[Boolean] = validate {
      backendConfig.as[Option[Boolean]]("batch.restrict-metadata-access").getOrElse(false)
    }
    val batchEnableFuse: ErrorOr[Boolean] = validate {
      backendConfig.as[Option[Boolean]]("batch.enable-fuse").getOrElse(false)
    }
    val logsPolicy: ErrorOr[GcpBatchLogsPolicy] = validate {
      backendConfig.as[Option[String]]("batch.logs-policy").getOrElse("CLOUD_LOGGING") match {
        case "CLOUD_LOGGING" => GcpBatchLogsPolicy.CloudLogging
        case "PATH" => GcpBatchLogsPolicy.Path
        case other =>
          throw new IllegalArgumentException(
            s"Unrecognized logs policy entry: $other. Supported strategies are CLOUD_LOGGING and PATH."
          )
      }
    }

    val dockerhubToken: ErrorOr[String] = validate {
      backendConfig.as[Option[String]]("dockerhub.token").getOrElse("")
    }

    val gcsFilesystemAuthName: ErrorOr[String] = validate {
      backendConfig.as[String]("filesystems.gcs.auth")
    }
    val qpsValidation = validateQps(backendConfig)
    val duplicationStrategy = validate {
      backendConfig.as[Option[String]]("filesystems.gcs.caching.duplication-strategy").getOrElse("copy") match {
        case "copy" => CopyCachedOutputs
        case "reference" => UseOriginalCachedOutputs
        case other =>
          throw new IllegalArgumentException(
            s"Unrecognized caching duplication strategy: $other. Supported strategies are copy and reference. See reference.conf for more details."
          )
      }
    }

    val requestWorkers: ErrorOr[Int Refined Positive] =
      validatePositiveInt(backendConfig.as[Option[Int]]("request-workers").getOrElse(3), "request-workers")

    val batchTimeout: FiniteDuration = backendConfig.getOrElse("batch-timeout", 7.days)

    val dockerMirroring: Option[DockerMirroring] = DockerMirroring.fromConfig(backendConfig)

    val logFlushPeriod: Option[FiniteDuration] = backendConfig.as[Option[FiniteDuration]]("log-flush-period") match {
      case Some(duration) if duration.isFinite => Option(duration)
      // "Inf" disables upload
      case Some(_) => None
      // Defaults to 1 minute
      case None => Option(1.minute)
    }

    val parallelCompositeUploadThreshold =
      validateGsutilMemorySpecification(backendConfig, "batch.parallel-composite-upload-threshold")

    val localizationAttempts: ErrorOr[Int Refined Positive] = backendConfig
      .as[Option[Int]]("batch.localization-attempts")
      .map(attempts => validatePositiveInt(attempts, "batch.localization-attempts"))
      .getOrElse(DefaultGcsTransferAttempts.validNel)

    val gcsTransferConfiguration: ErrorOr[GcsTransferConfiguration] =
      (localizationAttempts, parallelCompositeUploadThreshold) mapN GcsTransferConfiguration.apply

    val vpcNetworkName: ErrorOr[Option[String]] = validate {
      backendConfig.getAs[String]("virtual-private-cloud.network-name")
    }
    val vpcSubnetworkName: ErrorOr[Option[String]] = validate {
      backendConfig.getAs[String]("virtual-private-cloud.subnetwork-name")
    }
    val vpcNetworkLabel: ErrorOr[Option[String]] = validate {
      backendConfig.getAs[String]("virtual-private-cloud.network-label-key")
    }
    val vpcSubnetworkLabel: ErrorOr[Option[String]] = validate {
      backendConfig.getAs[String]("virtual-private-cloud.subnetwork-label-key")
    }
    val vpcAuth: ErrorOr[Option[String]] = validate {
      backendConfig.getAs[String]("virtual-private-cloud.auth")
    }

    val virtualPrivateCloudConfiguration: ErrorOr[VirtualPrivateCloudConfiguration] =
      (vpcNetworkName, vpcSubnetworkName, vpcNetworkLabel, vpcSubnetworkLabel, vpcAuth) flatMapN validateVPCConfig

    val batchRequestsReadTimeout =
      readOptionalPositiveMillisecondsIntFromDuration(backendConfig, "batch-requests.timeouts.read")
    val batchRequestsConnectTimeout =
      readOptionalPositiveMillisecondsIntFromDuration(backendConfig, "batch-requests.timeouts.connect")

    val batchRequestTimeoutConfigurationValidation = (batchRequestsReadTimeout, batchRequestsConnectTimeout) mapN {
      (read, connect) =>
        BatchRequestTimeoutConfiguration(readTimeoutMillis = read, connectTimeoutMillis = connect)
    }

    val referenceDiskLocalizationManifestFiles: ErrorOr[Option[List[ManifestFile]]] =
      validateReferenceDiskManifestConfigs(backendConfig, backendName)

    val referenceDiskValidationAuthNameOption: ErrorOr[Option[String]] = validate {
      backendConfig.getAs[String]("reference-disk-localization.validation.auth")
    }

    val checkpointingInterval: FiniteDuration = backendConfig.getOrElse(checkpointingIntervalKey, 10.minutes)

    val maxTransientErrorRetries: Int =
      backendConfig.as[Option[Int]]("max-transient-error-retries").getOrElse(10)

    def authGoogleConfigForBatchConfigurationAttributes(
      project: String,
      bucket: String,
      batchName: String,
      location: String,
      restrictMetadata: Boolean,
      dockerhubToken: String,
      enableFuse: Boolean,
      gcsName: String,
      qps: Int Refined Positive,
      cacheHitDuplicationStrategy: BatchCacheHitDuplicationStrategy,
      requestWorkers: Int Refined Positive,
      gcsTransferConfiguration: GcsTransferConfiguration,
      virtualPrivateCloudConfiguration: VirtualPrivateCloudConfiguration,
      batchRequestTimeoutConfiguration: BatchRequestTimeoutConfiguration,
      referenceDiskLocalizationManifestFilesOpt: Option[List[ManifestFile]],
      referenceDiskValidationAuthNameOption: Option[String],
      logsPolicy: GcpBatchLogsPolicy
    ): ErrorOr[GcpBatchConfigurationAttributes] = {
      val referenceDiskValidationAuthName = referenceDiskValidationAuthNameOption.getOrElse(batchName)
      (
        googleConfig.auth(batchName),
        googleConfig.auth(gcsName),
        googleConfig.auth(referenceDiskValidationAuthName)
      ) mapN { (batchAuth, gcsAuth, referenceDiskValidationAuth) =>
        val generatedReferenceFilesMappingOpt = referenceDiskLocalizationManifestFilesOpt map {
          generateReferenceFilesMapping(referenceDiskValidationAuth, _)
        }
        models.GcpBatchConfigurationAttributes(
          project = project,
          computeServiceAccount = computeServiceAccount,
          auths = GcpBatchAuths(batchAuth, gcsAuth),
          dockerhubToken = dockerhubToken,
          enableFuse = enableFuse,
          executionBucket = bucket,
          location = location,
          maxPollingInterval = maxPollingInterval,
          qps = qps,
          cacheHitDuplicationStrategy = cacheHitDuplicationStrategy,
          requestWorkers = requestWorkers,
          batchTimeout = batchTimeout,
          dockerMirroringOpt = dockerMirroring,
          logFlushPeriod = logFlushPeriod,
          gcsTransferConfiguration = gcsTransferConfiguration,
          virtualPrivateCloudConfiguration = virtualPrivateCloudConfiguration,
          batchRequestTimeoutConfiguration = batchRequestTimeoutConfiguration,
          referenceFileToDiskImageMappingOpt = generatedReferenceFilesMappingOpt,
          checkpointingInterval = checkpointingInterval,
          logsPolicy = logsPolicy,
          maxTransientErrorRetries = maxTransientErrorRetries
        )
      }
    }

    (project,
     executionBucket,
     batchAuthName,
     location,
     batchRestrictMetadataAccess,
     dockerhubToken,
     batchEnableFuse,
     gcsFilesystemAuthName,
     qpsValidation,
     duplicationStrategy,
     requestWorkers,
     gcsTransferConfiguration,
     virtualPrivateCloudConfiguration,
     batchRequestTimeoutConfigurationValidation,
     referenceDiskLocalizationManifestFiles,
     referenceDiskValidationAuthNameOption,
     logsPolicy
    ) flatMapN authGoogleConfigForBatchConfigurationAttributes match {
      case Valid(r) => r
      case Invalid(f) =>
        throw new IllegalArgumentException with MessageAggregation {
          override val exceptionContext = "Google Cloud Batch configuration is not valid: Errors"
          override val errorMessages: List[String] = f.toList
        }
    }
  }

  /**
    * Validate that the entries corresponding to "reference-disk-localization-manifests" in the specified
    * backend are parseable as `ManifestFile`s.
    */
  private[batch] def validateReferenceDiskManifestConfigs(backendConfig: Config,
                                                          backendName: String
  ): ErrorOr[Option[List[ManifestFile]]] =
    Try(backendConfig.getAs[List[Config]]("reference-disk-localization-manifests")) match {
      case Failure(e) =>
        ("Error attempting to parse value for 'reference-disk-localization-manifests' as List[Config]: " +
          e.getMessage).invalidNel
      case Success(s) =>
        s match {
          case Some(configs) =>
            import _root_.io.circe.config.parser
            import _root_.io.circe.generic.auto._

            // Unfortunately the `as` method of `config` collides with the Ficus method of the same name, so invoke its
            // equivalent using clunkier syntax:
            configs traverse parser.decode[ManifestFile] match {
              case Right(manifests) =>
                logger.info(
                  s"Reference disks feature for $backendName backend is configured with the following reference images: ${manifests.map(_.imageIdentifier).mkString(", ")}."
                )
                Option(manifests).validNel
              case Left(err) =>
                val message =
                  s"Reference disks misconfigured for backend $backendName, could not parse as List[ManifestFile]"
                logger.error(message, err.getCause)
                s"$message: ${err.getMessage}".invalidNel
            }
          case None =>
            logger.info(s"Reference disks feature for $backendName backend is not configured.")
            None.validNel
        }
    }

  def validateQps(config: Config): ErrorOr[Int Refined Positive] = {
    import eu.timepit.refined._

    val qp100s = config.as[Option[Int]]("batch-queries-per-100-seconds").getOrElse(BatchApiDefaultQp100s)
    val qpsCandidate = qp100s / 100

    refineV[Positive](qpsCandidate) match {
      case Left(_) =>
        s"Calculated QPS for Google Batch API ($qpsCandidate/s) was not a positive integer (supplied value was $qp100s per 100s)".invalidNel
      case Right(refined) => refined.validNel
    }
  }

  def validateGsutilMemorySpecification(config: Config, configPath: String): ErrorOr[String] = {
    val entry = config.as[Option[String]](configPath)
    entry match {
      case None => "0".validNel
      case Some(v @ GsutilHumanBytes(_, _)) => v.validNel
      case Some(bad) =>
        s"Invalid gsutil memory specification in Cromwell configuration at path '$configPath': '$bad'".invalidNel
    }
  }

  def validatePositiveInt(n: Int, configPath: String): Validated[NonEmptyList[String], Refined[Int, Positive]] =
    refineV[Positive](n) match {
      case Left(_) => s"Value $n for $configPath is not strictly positive".invalidNel
      case Right(refined) => refined.validNel
    }

  def readOptionalPositiveMillisecondsIntFromDuration(backendConfig: Config,
                                                      configPath: String
  ): ErrorOr[Option[Int Refined Positive]] = {

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

  // Copy/port of gsutil's "_GenerateSuffixRegex"
  private[batch] lazy val GsutilHumanBytes: Regex = {
    val _EXP_STRINGS = List(
      List("B", "bit"),
      List("KiB", "Kibit", "K"),
      List("MiB", "Mibit", "M"),
      List("GiB", "Gibit", "G"),
      List("TiB", "Tibit", "T"),
      List("PiB", "Pibit", "P"),
      List("EiB", "Eibit", "E")
    )

    val suffixes = for {
      unit <- _EXP_STRINGS
      name <- unit
    } yield name

    // Differs from the Python original in a couple of ways:
    // * The Python original uses named groups which are not supported in Scala regexes.
    //   (?P<num>\d*\.\d+|\d+)\s*(?P<suffix>%s)?
    //
    // * The Python original lowercases both the units and the human string before running the matcher.
    //   This Scala version turns on the (?i) case insensitive matching regex option instead.
    val orSuffixes = suffixes.mkString("|")
    "(?i)(\\d*\\.\\d+|\\d+)\\s*(%s)?".format(orSuffixes).r
  }
}
