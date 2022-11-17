package cromwell.backend.impl.tes

import cats.syntax.validated._
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import wom.RuntimeAttributesKeys
import wom.format.MemorySize
import wom.types.{WomIntegerType, WomStringType}
import wom.values._

case class TesRuntimeAttributes(continueOnReturnCode: ContinueOnReturnCode,
                                dockerImage: String,
                                dockerWorkingDir: Option[String],
                                failOnStderr: Boolean,
                                cpu: Option[Int Refined Positive],
                                memory: Option[MemorySize],
                                disk: Option[MemorySize],
                                preemptible: Boolean,
                                backendParameters: Map[String, Option[String]])

object TesRuntimeAttributes {

  val DockerWorkingDirKey = "dockerWorkingDir"
  val DiskSizeKey = "disk"
  val PreemptibleKey = "preemptible"

  private def cpuValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Int Refined Positive] = CpuValidation.optional

  private def failOnStderrValidation(runtimeConfig: Option[Config]) = FailOnStderrValidation.default(runtimeConfig)

  private def continueOnReturnCodeValidation(runtimeConfig: Option[Config]) = ContinueOnReturnCodeValidation.default(runtimeConfig)

  private def diskSizeValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[MemorySize] = MemoryValidation.optional(DiskSizeKey)

  private def memoryValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[MemorySize] = MemoryValidation.optional(RuntimeAttributesKeys.MemoryKey)

  private val dockerValidation: RuntimeAttributesValidation[String] = DockerValidation.instance

  private val dockerWorkingDirValidation: OptionalRuntimeAttributesValidation[String] = DockerWorkingDirValidation.optional

  private def preemptibleValidation(runtimeConfig: Option[Config]) = PreemptibleValidation.default(runtimeConfig)

  def runtimeAttributesBuilder(backendRuntimeConfig: Option[Config]): StandardValidatedRuntimeAttributesBuilder =
    StandardValidatedRuntimeAttributesBuilder.default(backendRuntimeConfig).withValidation(
      cpuValidation(backendRuntimeConfig),
      memoryValidation(backendRuntimeConfig),
      diskSizeValidation(backendRuntimeConfig),
      dockerValidation,
      dockerWorkingDirValidation,
      preemptibleValidation(backendRuntimeConfig),
    )

  def makeBackendParameters(runtimeAttributes: Map[String, WomValue],
                            keysToExclude: Set[String],
                            config: TesConfiguration): Map[String, Option[String]] = {

    if (config.useBackendParameters)
      runtimeAttributes
        .view.filterKeys(k => !keysToExclude.contains(k))
        .flatMap( _ match {
          case (key, WomString(s)) => Option((key, Option(s)))
          case (key, WomOptionalValue(WomStringType, Some(WomString(optS)))) => Option((key, Option(optS)))
          case (key, WomOptionalValue(WomStringType, None)) => Option((key, None))
          case _ => None
        }).toMap
    else
      Map.empty
  }

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes, rawRuntimeAttributes: Map[String, WomValue], config: TesConfiguration): TesRuntimeAttributes = {
    val backendRuntimeConfig = config.runtimeConfig
    val docker: String = RuntimeAttributesValidation.extract(dockerValidation, validatedRuntimeAttributes)
    val dockerWorkingDir: Option[String] = RuntimeAttributesValidation.extractOption(dockerWorkingDirValidation.key, validatedRuntimeAttributes)
    val cpu: Option[Int Refined Positive] = RuntimeAttributesValidation.extractOption(cpuValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val memory: Option[MemorySize] = RuntimeAttributesValidation.extractOption(memoryValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val disk: Option[MemorySize] = RuntimeAttributesValidation.extractOption(diskSizeValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val failOnStderr: Boolean =
      RuntimeAttributesValidation.extract(failOnStderrValidation(backendRuntimeConfig), validatedRuntimeAttributes)
    val continueOnReturnCode: ContinueOnReturnCode =
      RuntimeAttributesValidation.extract(continueOnReturnCodeValidation(backendRuntimeConfig), validatedRuntimeAttributes)
    val preemptible: Boolean =
      RuntimeAttributesValidation.extract(preemptibleValidation(backendRuntimeConfig), validatedRuntimeAttributes)

    // !! NOTE !! If new validated attributes are added to TesRuntimeAttributes, be sure to include
    // their validations here so that they will be handled correctly with backendParameters.
    val validations = Set(
      dockerValidation,
      dockerWorkingDirValidation,
      cpuValidation(backendRuntimeConfig),
      memoryValidation(backendRuntimeConfig),
      diskSizeValidation(backendRuntimeConfig),
      failOnStderrValidation(backendRuntimeConfig),
      continueOnReturnCodeValidation(backendRuntimeConfig),
      preemptibleValidation(backendRuntimeConfig)
    )

    // BT-458 any strings included in runtime attributes that aren't otherwise used should be
    // passed through to the TES server as part of backend_parameters
    val keysToExclude = validations map { _.key }
    val backendParameters = makeBackendParameters(rawRuntimeAttributes, keysToExclude, config)

    new TesRuntimeAttributes(
      continueOnReturnCode,
      docker,
      dockerWorkingDir,
      failOnStderr,
      cpu,
      memory,
      disk,
      preemptible,
      backendParameters
    )
  }
}

object DockerWorkingDirValidation {
  lazy val instance: RuntimeAttributesValidation[String] = new DockerWorkingDirValidation
  lazy val optional: OptionalRuntimeAttributesValidation[String] = instance.optional
}

class DockerWorkingDirValidation extends StringRuntimeAttributesValidation(TesRuntimeAttributes.DockerWorkingDirKey) {
  // NOTE: Docker's current test specs don't like WdlInteger, etc. auto converted to WdlString.
  override protected def validateValue: PartialFunction[WomValue, ErrorOr[String]] = {
    case WomString(value) => value.validNel
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
  def default(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Boolean] = instance.withDefault(
    configDefaultWdlValue(runtimeConfig) getOrElse WomBoolean(false))
  def configDefaultWdlValue(runtimeConfig: Option[Config]): Option[WomValue] = instance.configDefaultWomValue(runtimeConfig)
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
