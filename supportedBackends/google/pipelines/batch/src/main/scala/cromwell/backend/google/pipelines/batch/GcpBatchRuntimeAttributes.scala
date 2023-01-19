package cromwell.backend.google.pipelines.batch

import com.typesafe.config.Config
import cromwell.backend.validation.BooleanRuntimeAttributesValidation
import wom.values.WomBoolean
import cromwell.backend.validation.{DockerValidation, StringRuntimeAttributesValidation, ValidatedRuntimeAttributes}
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation.{OptionalRuntimeAttributesValidation, RuntimeAttributesValidation}
import eu.timepit.refined.api.Refined
//import eu.timepit.refined.numeric
import cromwell.backend.validation._
import eu.timepit.refined.numeric.Positive


case class GcpBatchRuntimeAttributes(
                                      cpu: Int Refined Positive,
                                      cpuPlatform: Option[String],
                                      dockerImage: String,
                                      noAddress: Boolean)

object GcpBatchRuntimeAttributes {

  val NoAddressKey = "noAddress"
  private val noAddressValidationInstance = new BooleanRuntimeAttributesValidation(NoAddressKey)
  private val NoAddressDefaultValue = WomBoolean(false)


  val CpuPlatformKey = "cpuPlatform"
  private val cpuPlatformValidationInstance = new StringRuntimeAttributesValidation(CpuPlatformKey)
    .optional
  // via `gcloud compute zones describe us-central1-a`
  //val CpuPlatformIntelCascadeLakeValue = "Intel Cascade Lake"
  //val CpuPlatformAMDRomeValue = "AMD Rome"

  private def cpuValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int Refined Positive] = CpuValidation
    .instance
    .withDefault(CpuValidation
      .configDefaultWomValue(runtimeConfig) getOrElse CpuValidation
      .defaultMin)
  private def cpuPlatformValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] = cpuPlatformValidationInstance
  private val dockerValidation: RuntimeAttributesValidation[String] = DockerValidation.instance

  private def noAddressValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Boolean] = noAddressValidationInstance
    .withDefault(noAddressValidationInstance
      .configDefaultWomValue(runtimeConfig) getOrElse NoAddressDefaultValue)

  def runtimeAttributesBuilder(batchConfiguration: GcpBatchConfiguration): StandardValidatedRuntimeAttributesBuilder = {
    val runtimeConfig = batchConfiguration.runtimeConfig
    StandardValidatedRuntimeAttributesBuilder.default(runtimeConfig).withValidation(
      cpuValidation(runtimeConfig),
      cpuPlatformValidation(runtimeConfig),
      noAddressValidation(runtimeConfig),
      dockerValidation
    )
  }

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes, runtimeAttrsConfig: Option[Config]): GcpBatchRuntimeAttributes = {
    val cpu: Int Refined Positive = RuntimeAttributesValidation.extract(cpuValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val cpuPlatform: Option[String] = RuntimeAttributesValidation.extractOption(cpuPlatformValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)
    val docker: String = RuntimeAttributesValidation.extract(dockerValidation, validatedRuntimeAttributes)
    val noAddress: Boolean = RuntimeAttributesValidation.extract(noAddressValidation(runtimeAttrsConfig), validatedRuntimeAttributes)

    new GcpBatchRuntimeAttributes(
      cpu,
      cpuPlatform,
      docker,
      noAddress)
  }


}
