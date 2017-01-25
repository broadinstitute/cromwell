package cromwell.backend.impl.tes

import cats.syntax.validated._
import cromwell.backend.MemorySize
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation._
import lenthall.validation.ErrorOr.ErrorOr
import org.slf4j.Logger
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
}

class DiskSizeValidation extends MemoryValidation {
  override def key = TesRuntimeAttributes.DiskSizeKey
}


