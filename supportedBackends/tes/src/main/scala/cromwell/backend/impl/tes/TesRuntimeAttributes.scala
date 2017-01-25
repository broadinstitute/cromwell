package cromwell.backend.impl.tes

import cats.syntax.validated._
import cromwell.backend.MemorySize
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation._
import lenthall.validation.ErrorOr.ErrorOr
import org.slf4j.Logger
import wdl4s.parser.MemoryUnit
import wdl4s.values.{WdlInteger, WdlString, WdlValue}

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

  private val cpuValidation: RuntimeAttributesValidation[Int] = CpuValidation.default

  private val diskSizeValidation: RuntimeAttributesValidation[MemorySize] =
    DiskSizeValidation.withDefaultDiskSize(MemorySize.parse(DiskSizeDefaultValue).get)

  private val dockerValidation: RuntimeAttributesValidation[String] = DockerValidation.instance

  private val dockerWorkingDirValidation: OptionalRuntimeAttributesValidation[String] = DockerWorkingDirValidation.optional

  private val memoryValidation: RuntimeAttributesValidation[MemorySize] =
    MemoryValidation.withDefaultMemory(MemorySize.parse(MemoryDefaultValue).get)

  def runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    StandardValidatedRuntimeAttributesBuilder.default.withValidation(
      cpuValidation,
      memoryValidation,
      diskSizeValidation,
      dockerValidation,
      dockerWorkingDirValidation
    )

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes): TesRuntimeAttributes = {
    val docker: String = RuntimeAttributesValidation.extract(dockerValidation, validatedRuntimeAttributes)
    val dockerWorkingDir: Option[String] = RuntimeAttributesValidation.extractOption(dockerWorkingDirValidation.key, validatedRuntimeAttributes)
    val cpu: Int = RuntimeAttributesValidation.extract(cpuValidation, validatedRuntimeAttributes)
    val memory: MemorySize = RuntimeAttributesValidation.extract(memoryValidation, validatedRuntimeAttributes)
    val disk: MemorySize = RuntimeAttributesValidation.extract(diskSizeValidation, validatedRuntimeAttributes)
    val failOnStderr: Boolean =
      RuntimeAttributesValidation.extract(FailOnStderrValidation.default, validatedRuntimeAttributes)
    val continueOnReturnCode: ContinueOnReturnCode =
      RuntimeAttributesValidation.extract(ContinueOnReturnCodeValidation.default, validatedRuntimeAttributes)

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

  // NOTE: Currently only used by test specs
  private[tes] def apply(attrs: Map[String, WdlValue],
                         logger: Logger): TesRuntimeAttributes = {
    val runtimeAttributesBuilder = TesRuntimeAttributes.runtimeAttributesBuilder
    val validatedRuntimeAttributes = runtimeAttributesBuilder.build(attrs, logger)
    apply(validatedRuntimeAttributes)
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

  def withDefaultDiskSize(memorySize: MemorySize): RuntimeAttributesValidation[MemorySize] =
    instance.withDefault(WdlInteger(memorySize.bytes.toInt))

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


