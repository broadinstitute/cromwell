package cromwell.backend.impl.tes

import cats.syntax.validated._
import cromwell.backend.MemorySize
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation._
import lenthall.validation.ErrorOr.ErrorOr
import org.slf4j.Logger
import wdl4s.values.{WdlString, WdlValue}

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

  val DiskKey = "disk"
  private val DiskDefaultValue = "2 GB"

  private val cpuValidation: RuntimeAttributesValidation[Int] = CpuValidation.default

  private val diskSizeValidation: RuntimeAttributesValidation[MemorySize] =
    MemoryValidation.withDefaultMemory(MemorySize.parse(DiskDefaultValue).get)

  private val dockerValidation: RuntimeAttributesValidation[String] = DockerValidation.instance

  private val dockerWorkingDirValidation: RuntimeAttributesValidation[String] = DockerWorkingDirValidation.instance

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
    val dockerWorkingDir: Option[String] = RuntimeAttributesValidation.extractOption(dockerWorkingDirValidation, validatedRuntimeAttributes)
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
  override protected def missingValueMessage: String = s"Expecting $key runtime attribute to be an String"

  override protected def validateValue: PartialFunction[WdlValue, ErrorOr[String]] = {
    case WdlString(value) => value.validNel
  }
}
