package cromwell.backend.impl.tes

import cats.data.Validated
import cats.syntax.validated._
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.google.batch.io.{BatchApiEmptyMountedDisk, GcpBatchAttachedDisk, GcpBatchWorkingDisk}
import cromwell.backend.google.batch.models.DisksValidation
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import wdl4s.parser.MemoryUnit
import wom.RuntimeAttributesKeys
import wom.format.MemorySize
import wom.types.{WomIntegerType, WomStringType}
import wom.values._

import java.util.regex.Pattern
import org.slf4j.LoggerFactory

case class TesRuntimeAttributes(continueOnReturnCode: ContinueOnReturnCode,
                                dockerImage: String,
                                dockerWorkingDir: Option[String],
                                failOnStderr: Boolean,
                                cpu: Option[Int Refined Positive],
                                memory: Option[MemorySize],
                                disk: Option[MemorySize],
                                preemptible: Boolean,
                                localizedSasEnvVar: Option[String],
                                backendParameters: Map[String, Option[String]],
                                memoryRetryMultiplier: Option[Double]
)

object TesRuntimeAttributes {
  private val log = LoggerFactory.getLogger(getClass.getSimpleName)

  val DockerWorkingDirKey = "dockerWorkingDir"
  val DiskSizeKey = "disk"
  val PreemptibleKey = "preemptible"
  val LocalizedSasKey = "azureSasEnvironmentVariable"
  val MemoryRetryMultiplierKey = "memory_retry_multiplier"
  val BackoffLimitKey = "backoff_limit"

  private def cpuValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Int Refined Positive] =
    CpuValidation.configDefaultWomValue(runtimeConfig) match {
      case Some(default) => CpuValidation.instance.withDefault(default).optional
      case None          => CpuValidation.optional
    }

  private def failOnStderrValidation(runtimeConfig: Option[Config]) = FailOnStderrValidation.default(runtimeConfig)

  private def continueOnReturnCodeValidation(runtimeConfig: Option[Config]) =
    ContinueOnReturnCodeValidation.default(runtimeConfig)

  private def diskSizeValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[MemorySize] =
    MemoryValidation.configDefaultString(DiskSizeKey, runtimeConfig) match {
      case Some(default) => MemoryValidation.withDefaultMemory(DiskSizeKey, default).optional
      case None          => MemoryValidation.optional(DiskSizeKey)
    }

  private def diskSizeCompatValidation(
    runtimeConfig: Option[Config]
  ): OptionalRuntimeAttributesValidation[Seq[GcpBatchAttachedDisk]] =
    DisksValidation.optional

  private def memoryValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[MemorySize] =
    MemoryValidation.configDefaultString(RuntimeAttributesKeys.MemoryKey, runtimeConfig) match {
      case Some(default) => MemoryValidation.withDefaultMemory(RuntimeAttributesKeys.MemoryKey, default).optional
      case None          => MemoryValidation.optional(RuntimeAttributesKeys.MemoryKey)
    }

  // As of WDL 1.1 these two are aliases of each other
  private val dockerValidation: OptionalRuntimeAttributesValidation[Containers] = DockerValidation.instance
  private val containerValidation: OptionalRuntimeAttributesValidation[Containers] = ContainerValidation.instance

  private val dockerWorkingDirValidation: OptionalRuntimeAttributesValidation[String] =
    DockerWorkingDirValidation.optional
  private def preemptibleValidation(runtimeConfig: Option[Config]) = PreemptibleValidation.default(runtimeConfig)
  private def localizedSasValidation: OptionalRuntimeAttributesValidation[String] = LocalizedSasValidation.optional

  private def memoryRetryMultiplierValidation(
    runtimeConfig: Option[Config]
  ): OptionalRuntimeAttributesValidation[Double] = {
    val instance = new FloatRuntimeAttributesValidation(MemoryRetryMultiplierKey)
    instance.configDefaultWomValue(runtimeConfig) match {
      case Some(default) => instance.withDefault(default).optional
      case None          => instance.optional
    }
  }

  private def backoffLimitValidation(
    runtimeConfig: Option[Config]
  ): OptionalRuntimeAttributesValidation[String] = {
    val instance = new BackoffLimitValidation
    instance.configDefaultWomValue(runtimeConfig) match {
      case Some(default) => instance.withDefault(default).optional
      case None          => instance.optional
    }
  }

  def runtimeAttributesBuilder(backendRuntimeConfig: Option[Config]): StandardValidatedRuntimeAttributesBuilder =
    // !! NOTE !! If new validated attributes are added to TesRuntimeAttributes, be sure to include
    // their validations here so that they will be handled correctly with backendParameters.
    // Location 2 of 2
    StandardValidatedRuntimeAttributesBuilder
      .default(backendRuntimeConfig)
      .withValidation(
        cpuValidation(backendRuntimeConfig),
        memoryValidation(backendRuntimeConfig),
        diskSizeValidation(backendRuntimeConfig),
        diskSizeCompatValidation(backendRuntimeConfig),
        dockerValidation,
        containerValidation,
        dockerWorkingDirValidation,
        preemptibleValidation(backendRuntimeConfig),
        localizedSasValidation,
        memoryRetryMultiplierValidation(backendRuntimeConfig),
        backoffLimitValidation(backendRuntimeConfig)
      )

  def makeBackendParameters(runtimeAttributes: Map[String, WomValue],
                            keysToExclude: Set[String],
                            config: TesConfiguration
  ): Map[String, Option[String]] = {
    // unknownKeys was declared but never used (compiler error with -Wunused).
    // Log at debug so callers can see which runtime attributes are being forwarded as
    // TES backend_parameters (useful for diagnosing Funnel/TES 1.1 passthrough issues).
    val unknownKeys = runtimeAttributes.keySet -- keysToExclude
    if (unknownKeys.nonEmpty)
      log.debug("makeBackendParameters: forwarding non-standard runtime keys as TES backend_parameters: {}",
                unknownKeys.mkString(", "))

    if (config.useBackendParameters)
      runtimeAttributes.view
        .filterKeys(k => !keysToExclude.contains(k))
        .flatMap(_ match {
          case (key, WomString(s)) => Option((key, Option(s)))
          case (key, WomOptionalValue(WomStringType, Some(WomString(optS)))) => Option((key, Option(optS)))
          case (key, WomOptionalValue(WomStringType, None)) => Option((key, None))
          case _ => None
        })
        .toMap
    else
      Map.empty
  }

  private def detectDiskFormat(backendRuntimeConfig: Option[Config],
                               validatedRuntimeAttributes: ValidatedRuntimeAttributes
  ): Option[MemorySize] = {

    def adaptPapiDisks(disks: Seq[GcpBatchAttachedDisk]): MemorySize =
      disks match {
        case disk :: Nil if disk.isInstanceOf[GcpBatchWorkingDisk] =>
          MemorySize(disk.sizeGb.toDouble, MemoryUnit.GB)
        case _ :: _ =>
          // When a user specifies only a custom disk, we add the default disk in the background, so we technically have multiple disks.
          // But we don't want to confuse the user with `multiple disks` message when they only put one.
          if (disks.exists(_.isInstanceOf[BatchApiEmptyMountedDisk]))
            throw new IllegalArgumentException("Disks with custom mount points are not supported by this backend")
          else
            // Multiple `local-disk` is not legal, but possible and should be detected
            throw new IllegalArgumentException("Expecting exactly one disk definition on this backend, found multiple")
      }

    val maybeTesDisk: Option[MemorySize] =
      RuntimeAttributesValidation.extractOption(diskSizeValidation(backendRuntimeConfig).key,
                                                validatedRuntimeAttributes
      )
    val maybePapiDisk: Option[Seq[GcpBatchAttachedDisk]] =
      RuntimeAttributesValidation.extractOption(diskSizeCompatValidation(backendRuntimeConfig).key,
                                                validatedRuntimeAttributes
      )

    (maybeTesDisk, maybePapiDisk) match {
      case (Some(tesDisk: MemorySize), _) =>
        Option(
          tesDisk
        ) // If WDLs are in circulation with both `disk` and `disks`, pick the one intended for this backend
      case (None, Some(papiDisks: Seq[GcpBatchAttachedDisk])) =>
        Option(adaptPapiDisks(papiDisks))
      case _ =>
        None
    }
  }

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes,
            rawRuntimeAttributes: Map[String, WomValue],
            config: TesConfiguration
  ): TesRuntimeAttributes = {
    val backendRuntimeConfig = config.runtimeConfig
    val docker: String = Containers.extractContainer(validatedRuntimeAttributes)
    val dockerWorkingDir: Option[String] =
      RuntimeAttributesValidation.extractOption(dockerWorkingDirValidation.key, validatedRuntimeAttributes)
    val cpu: Option[Int Refined Positive] =
      RuntimeAttributesValidation.extractOption(cpuValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val memory: Option[MemorySize] =
      RuntimeAttributesValidation.extractOption(memoryValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val disk: Option[MemorySize] = detectDiskFormat(backendRuntimeConfig, validatedRuntimeAttributes)
    val failOnStderr: Boolean =
      RuntimeAttributesValidation.extract(failOnStderrValidation(backendRuntimeConfig), validatedRuntimeAttributes)
    val continueOnReturnCode: ContinueOnReturnCode =
      RuntimeAttributesValidation.extract(continueOnReturnCodeValidation(backendRuntimeConfig),
                                          validatedRuntimeAttributes
      )
    val preemptible: Boolean =
      RuntimeAttributesValidation.extract(preemptibleValidation(backendRuntimeConfig), validatedRuntimeAttributes)
    val localizedSas: Option[String] =
      RuntimeAttributesValidation.extractOption(localizedSasValidation.key, validatedRuntimeAttributes)

    // !! NOTE !! If new validated attributes are added to TesRuntimeAttributes, be sure to include
    // their validations here so that they will be handled correctly with backendParameters.
    // Location 1 of 2
    val validations = Set(
      dockerValidation,
      containerValidation,
      dockerWorkingDirValidation,
      cpuValidation(backendRuntimeConfig),
      memoryValidation(backendRuntimeConfig),
      diskSizeValidation(backendRuntimeConfig),
      diskSizeCompatValidation(backendRuntimeConfig),
      failOnStderrValidation(backendRuntimeConfig),
      continueOnReturnCodeValidation(backendRuntimeConfig),
      preemptibleValidation(backendRuntimeConfig),
      localizedSasValidation,
      memoryRetryMultiplierValidation(backendRuntimeConfig),
      backoffLimitValidation(backendRuntimeConfig)
    )

    val memoryRetryMultiplier: Option[Double] =
      RuntimeAttributesValidation.extractOption(memoryRetryMultiplierValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)

    // Extract backoff_limit from validated attributes (covers both default-runtime-attributes
    // config defaults and per-task runtime { backoff_limit: "N" } declarations).
    val backoffLimit: Option[String] =
      RuntimeAttributesValidation.extractOption(backoffLimitValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    println(s"[backoff_limit] apply: extracted backoffLimit from validatedRuntimeAttributes = $backoffLimit")
    println(s"[backoff_limit] apply: rawRuntimeAttributes keys = ${rawRuntimeAttributes.keys.mkString(", ")}")

    // BT-458 any strings included in runtime attributes that aren't otherwise used should be
    // passed through to the TES server as part of backend_parameters
    val keysToExclude = validations map { _.key }
    val rawBackendParameters = makeBackendParameters(rawRuntimeAttributes, keysToExclude, config)
    // Inject backoff_limit from validated attributes (handles the config-default path, since
    // default-runtime-attributes values don't appear in rawRuntimeAttributes for unknown keys).
    val backendParameters = if (config.useBackendParameters)
      backoffLimit.fold(rawBackendParameters)(v => rawBackendParameters + (BackoffLimitKey -> Option(v)))
    else rawBackendParameters
    println(s"[backoff_limit] apply: final backendParameters = $backendParameters")

    new TesRuntimeAttributes(
      continueOnReturnCode,
      docker,
      dockerWorkingDir,
      failOnStderr,
      cpu,
      memory,
      disk,
      preemptible,
      localizedSas,
      backendParameters,
      memoryRetryMultiplier
    )
  }
}

object DockerWorkingDirValidation {
  lazy val instance: RuntimeAttributesValidation[String] = new DockerWorkingDirValidation
  lazy val optional: OptionalRuntimeAttributesValidation[String] = instance.optional
}

class DockerWorkingDirValidation extends StringRuntimeAttributesValidation(TesRuntimeAttributes.DockerWorkingDirKey) {
  // NOTE: Docker's current test specs don't like WdlInteger, etc. auto converted to WdlString.
  override protected def validateValue: PartialFunction[WomValue, ErrorOr[String]] = { case WomString(value) =>
    value.validNel
  }
}

/**
  * Validates the "preemptible" runtime attribute as a Boolean or a String 'true' or 'false', returning the value as a
  * `Boolean`.
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * `configDefaultWdlValue` returns the value of the attribute as specified by the
  * reference.conf file, coerced into a WomValue.
  *
  * `default` a validation with the default value specified by the reference.conf file.
  */

object PreemptibleValidation {
  lazy val instance: RuntimeAttributesValidation[Boolean] = new PreemptibleValidation
  def default(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Boolean] =
    instance.withDefault(configDefaultWdlValue(runtimeConfig) getOrElse WomBoolean(false))
  def configDefaultWdlValue(runtimeConfig: Option[Config]): Option[WomValue] =
    instance.configDefaultWomValue(runtimeConfig)
}

class PreemptibleValidation extends BooleanRuntimeAttributesValidation(TesRuntimeAttributes.PreemptibleKey) {
  override def usedInCallCaching: Boolean = false

  override protected def validateExpression: PartialFunction[WomValue, Boolean] = {
    case womBoolValue if womType.coerceRawValue(womBoolValue).isSuccess => true
    case womIntValue if WomIntegerType.coerceRawValue(womIntValue).isSuccess => true
  }

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Boolean]] = {
    case value if womType.coerceRawValue(value).isSuccess =>
      validateCoercedValue(womType.coerceRawValue(value).get.asInstanceOf[WomBoolean])
    // The TES spec requires a boolean preemptible value, but many WDLs written originally
    // for other backends use an integer. Interpret integers > 0 as true, others as false.
    case value if WomIntegerType.coerceRawValue(value).isSuccess =>
      validateCoercedValue(WomBoolean(WomIntegerType.coerceRawValue(value).get.asInstanceOf[WomInteger].value > 0))
    case value if womType.coerceRawValue(value.valueString).isSuccess =>
      /*
      NOTE: This case statement handles WdlString("true") coercing to WdlBoolean(true).
      For some reason "true" as String is coercable... but not the WdlString.
       */
      validateCoercedValue(womType.coerceRawValue(value.valueString).get.asInstanceOf[WomBoolean])
  }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be an Integer, Boolean, or a String with values of 'true' or 'false'"
}

object BackoffLimitValidation {
  lazy val instance: RuntimeAttributesValidation[String] = new BackoffLimitValidation
  lazy val optional: OptionalRuntimeAttributesValidation[String] = instance.optional
}

class BackoffLimitValidation extends StringRuntimeAttributesValidation(TesRuntimeAttributes.BackoffLimitKey) {
  override protected def validateValue: PartialFunction[WomValue, ErrorOr[String]] = {
    case WomString(value)  => value.validNel
    case WomInteger(value) => value.toString.validNel
  }
}

object LocalizedSasValidation {
  lazy val instance: RuntimeAttributesValidation[String] = new LocalizedSasValidation
  lazy val optional: OptionalRuntimeAttributesValidation[String] = instance.optional
}

class LocalizedSasValidation extends StringRuntimeAttributesValidation(TesRuntimeAttributes.LocalizedSasKey) {
  private def isValidBashVariableName(str: String): Boolean = {
    // require string be only letters, numbers, and underscores
    val pattern = Pattern.compile("^[a-zA-Z0-9_]+$", Pattern.CASE_INSENSITIVE)
    val matcher = pattern.matcher(str)
    matcher.find
  }

  override protected def invalidValueMessage(value: WomValue): String =
    s"Invalid Runtime Attribute value for ${TesRuntimeAttributes.LocalizedSasKey}. Value must be a string containing only letters, numbers, and underscores."

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[String]] = { case WomString(value) =>
    if (isValidBashVariableName(value)) value.validNel else Validated.invalidNel(invalidValueMessage(WomString(value)))
  }
}
