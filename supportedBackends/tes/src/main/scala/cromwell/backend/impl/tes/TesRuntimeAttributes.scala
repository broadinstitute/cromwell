package cromwell.backend.impl.tes

import cats.syntax.validated._
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import squants.information.Information
import wom.RuntimeAttributesKeys
import wom.values._

case class TesRuntimeAttributes(continueOnReturnCode: ContinueOnReturnCode,
                                dockerImage: String,
                                dockerWorkingDir: Option[String],
                                failOnStderr: Boolean,
                                cpu: Option[Int Refined Positive],
                                memory: Option[Information],
                                disk: Option[Information])

object TesRuntimeAttributes {

  val DockerWorkingDirKey = "dockerWorkingDir"
  val DiskSizeKey = "disk"

  private def cpuValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Int Refined Positive] = CpuValidation.optional

  private def failOnStderrValidation(runtimeConfig: Option[Config]) = FailOnStderrValidation.default(runtimeConfig)

  private def continueOnReturnCodeValidation(runtimeConfig: Option[Config]) = ContinueOnReturnCodeValidation.default(runtimeConfig)

  private def diskSizeValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Information] = MemoryValidation.optional(DiskSizeKey)

  private def memoryValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Information] = MemoryValidation.optional(RuntimeAttributesKeys.MemoryKey)

  private val dockerValidation: RuntimeAttributesValidation[String] = DockerValidation.instance

  private val dockerWorkingDirValidation: OptionalRuntimeAttributesValidation[String] = DockerWorkingDirValidation.optional

  def runtimeAttributesBuilder(backendRuntimeConfig: Option[Config]): StandardValidatedRuntimeAttributesBuilder =
    StandardValidatedRuntimeAttributesBuilder.default(backendRuntimeConfig).withValidation(
      cpuValidation(backendRuntimeConfig),
      memoryValidation(backendRuntimeConfig),
      diskSizeValidation(backendRuntimeConfig),
      dockerValidation,
      dockerWorkingDirValidation
    )

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes, backendRuntimeConfig: Option[Config]): TesRuntimeAttributes = {
    val docker: String = RuntimeAttributesValidation.extract(dockerValidation, validatedRuntimeAttributes)
    val dockerWorkingDir: Option[String] = RuntimeAttributesValidation.extractOption(dockerWorkingDirValidation.key, validatedRuntimeAttributes)
    val cpu: Option[Int Refined Positive] = RuntimeAttributesValidation.extractOption(cpuValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val memory: Option[Information] = RuntimeAttributesValidation.extractOption(memoryValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val disk: Option[Information] = RuntimeAttributesValidation.extractOption(diskSizeValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
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
  override protected def validateValue: PartialFunction[WomValue, ErrorOr[String]] = {
    case WomString(value) => value.validNel
  }
}

