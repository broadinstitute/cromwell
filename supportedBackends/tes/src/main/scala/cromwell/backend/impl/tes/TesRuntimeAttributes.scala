package cromwell.backend.impl.tes

import cats.syntax.validated._
import com.typesafe.config.Config
import cromwell.backend.MemorySize
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation._
import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.parser.MemoryUnit
import wdl4s.values.{WdlInteger, WdlString, WdlValue}

import scala.util.{Failure, Success}

case class TesRuntimeAttributes(continueOnReturnCode: ContinueOnReturnCode,
                                dockerImage: String,
                                dockerWorkingDir: Option[String],
                                failOnStderr: Boolean,
                                cpu: Int,
                                memory: MemorySize,
                                disk: MemorySize)

object TesRuntimeAttributes {

  val DockerWorkingDirKey = "dockerWorkingDir"

  private val MemoryDefaultValue = "2 GB"

  val DiskSizeKey = "disk"
  private val DiskSizeDefaultValue = "2 GB"

  private def cpuValidation(runtimeConfig: Config): RuntimeAttributesValidation[Int] = CpuValidation.instance
    .withDefault(CpuValidation.configDefaultWdlValue(runtimeConfig) getOrElse CpuValidation.default)

  private def failOnStderrValidation(runtimeConfig: Config) = FailOnStderrValidation.default(runtimeConfig)

  private def continueOnReturnCodeValidation(runtimeConfig: Config) = ContinueOnReturnCodeValidation.default(runtimeConfig)

  private def diskSizeValidation(runtimeConfig: Config): RuntimeAttributesValidation[MemorySize] = DiskSizeValidation
    .withDefaultDiskSize(DiskSizeValidation.configDefaultString(runtimeConfig) getOrElse DiskSizeDefaultValue)

  private def memoryValidation(runtimeConfig: Config): RuntimeAttributesValidation[MemorySize] = {
    MemoryValidation.withDefaultMemory(MemoryValidation.configDefaultString(runtimeConfig) getOrElse MemoryDefaultValue)
  }

  private val dockerValidation: RuntimeAttributesValidation[String] = DockerValidation.instance

  private val dockerWorkingDirValidation: OptionalRuntimeAttributesValidation[String] = DockerWorkingDirValidation.optional

  def runtimeAttributesBuilder(backendRuntimeConfig: Config): StandardValidatedRuntimeAttributesBuilder =
    StandardValidatedRuntimeAttributesBuilder.default(backendRuntimeConfig).withValidation(
      cpuValidation(backendRuntimeConfig),
      memoryValidation(backendRuntimeConfig),
      diskSizeValidation(backendRuntimeConfig),
      dockerValidation,
      dockerWorkingDirValidation
    )

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes, backendRuntimeConfig: Config): TesRuntimeAttributes = {
    val docker: String = RuntimeAttributesValidation.extract(dockerValidation, validatedRuntimeAttributes)
    val dockerWorkingDir: Option[String] = RuntimeAttributesValidation.extractOption(dockerWorkingDirValidation.key, validatedRuntimeAttributes)
    val cpu: Int = RuntimeAttributesValidation.extract(cpuValidation(backendRuntimeConfig), validatedRuntimeAttributes)
    val memory: MemorySize = RuntimeAttributesValidation.extract(memoryValidation(backendRuntimeConfig), validatedRuntimeAttributes)
    val disk: MemorySize = RuntimeAttributesValidation.extract(diskSizeValidation(backendRuntimeConfig), validatedRuntimeAttributes)
    val failOnStderr: Boolean =
      RuntimeAttributesValidation.extract(failOnStderrValidation(backendRuntimeConfig), validatedRuntimeAttributes)
    val continueOnReturnCode: ContinueOnReturnCode =
      RuntimeAttributesValidation.extract(continueOnReturnCodeValidation(backendRuntimeConfig), validatedRuntimeAttributes)

    new TesRuntimeAttributes(
      continueOnReturnCode,
      docker,
      dockerWorkingDir,
      failOnStderr,
      cpu,
      memory,
      disk
    )
  }
}

object DockerWorkingDirValidation {
  lazy val instance: RuntimeAttributesValidation[String] = new DockerWorkingDirValidation
  lazy val optional: OptionalRuntimeAttributesValidation[String] = instance.optional
}

class DockerWorkingDirValidation extends StringRuntimeAttributesValidation(TesRuntimeAttributes.DockerWorkingDirKey) {
  // NOTE: Docker's current test specs don't like WdlInteger, etc. auto converted to WdlString.
  override protected def validateValue: PartialFunction[WdlValue, ErrorOr[String]] = {
    case WdlString(value) => value.validNel
  }
}

object DiskSizeValidation {
  lazy val instance: RuntimeAttributesValidation[MemorySize] = new DiskSizeValidation
  lazy val optional: OptionalRuntimeAttributesValidation[MemorySize] = instance.optional
  def configDefaultString(runtimeConfig: Config): Option[String] = instance.configDefaultValue(runtimeConfig)
  def withDefaultDiskSize(memorySize: String): RuntimeAttributesValidation[MemorySize] = {
    MemorySize.parse(memorySize) match {
      case Success(memory) => instance.withDefault(WdlInteger(memory.bytes.toInt))
      case Failure(_) => instance.withDefault(WdlString(memorySize.toString))
    }
  }

  private val wrongAmountFormat =
    s"Expecting ${TesRuntimeAttributes.DiskSizeKey} runtime attribute value greater than 0 but got %s"
  private val wrongTypeFormat =
    s"Expecting ${TesRuntimeAttributes.DiskSizeKey} runtime attribute to be an Integer or String with format '8 GB'." +
      s" Exception: %s"

  def validateDiskSizeString(wdlString: WdlString): ErrorOr[MemorySize] =
    validateDiskSizeString(wdlString.value)

  def validateDiskSizeString(value: String): ErrorOr[MemorySize] = {
    MemorySize.parse(value) match {
      case scala.util.Success(memorySize: MemorySize) if memorySize.amount > 0 =>
        memorySize.to(MemoryUnit.GB).validNel
      case scala.util.Success(memorySize: MemorySize) =>
        wrongAmountFormat.format(memorySize.amount).invalidNel
      case scala.util.Failure(throwable) =>
        wrongTypeFormat.format(throwable.getMessage).invalidNel
    }
  }

  def validateDiskSizeInteger(wdlInteger: WdlInteger): ErrorOr[MemorySize] =
    validateDiskSizeInteger(wdlInteger.value)

  def validateDiskSizeInteger(value: Int): ErrorOr[MemorySize] = {
    if (value <= 0)
      wrongAmountFormat.format(value).invalidNel
    else
      MemorySize(value.toDouble, MemoryUnit.Bytes).to(MemoryUnit.GB).validNel
  }
}

class DiskSizeValidation extends MemoryValidation {
  override def key = TesRuntimeAttributes.DiskSizeKey

  override protected def validateValue: PartialFunction[WdlValue, ErrorOr[MemorySize]] = {
    case WdlInteger(value) => DiskSizeValidation.validateDiskSizeInteger(value)
    case WdlString(value) => DiskSizeValidation.validateDiskSizeString(value)
  }

  override def missingValueMessage: String = DiskSizeValidation.wrongTypeFormat.format("Not supported WDL type value")
}


