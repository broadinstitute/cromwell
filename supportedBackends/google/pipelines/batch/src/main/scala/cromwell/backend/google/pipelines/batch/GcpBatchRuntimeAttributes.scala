package cromwell.backend.google.pipelines.batch

import com.typesafe.config.Config
import cromwell.backend.google.pipelines.common.ZonesValidation
import cromwell.backend.validation.BooleanRuntimeAttributesValidation
import wom.values.{WomBoolean, WomString}
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
                                      zones: Vector[String],
                                      dockerImage: String,
                                      noAddress: Boolean,
                                      failOnStderr: Boolean)

object GcpBatchRuntimeAttributes {

  val ZonesKey = "zones"
  private val ZonesDefaultValue = WomString("us-central1-b")


  val NoAddressKey = "noAddress"
  private val noAddressValidationInstance = new BooleanRuntimeAttributesValidation(NoAddressKey)
  private val NoAddressDefaultValue = WomBoolean(false)


  val CpuPlatformKey = "cpuPlatform"
  private val cpuPlatformValidationInstance = new StringRuntimeAttributesValidation(CpuPlatformKey)
    .optional
  // via `gcloud compute zones describe us-central1-a`
  val CpuPlatformIntelCascadeLakeValue = "Intel Cascade Lake"
  val CpuPlatformAMDRomeValue = "AMD Rome"

  private def cpuValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int Refined Positive] = CpuValidation
    .instance
    .withDefault(CpuValidation
      .configDefaultWomValue(runtimeConfig) getOrElse CpuValidation
      .defaultMin)
  private def cpuPlatformValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] = cpuPlatformValidationInstance
  private val dockerValidation: RuntimeAttributesValidation[String] = DockerValidation.instance

  private def failOnStderrValidation(runtimeConfig: Option[Config]) = FailOnStderrValidation.default(runtimeConfig)
  private def zonesValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Vector[String]] = ZonesValidation
    .withDefault(ZonesValidation
      .configDefaultWomValue(runtimeConfig) getOrElse ZonesDefaultValue)

  private def noAddressValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Boolean] = noAddressValidationInstance
    .withDefault(noAddressValidationInstance
      .configDefaultWomValue(runtimeConfig) getOrElse NoAddressDefaultValue)

  def runtimeAttributesBuilder(batchConfiguration: GcpBatchConfiguration): StandardValidatedRuntimeAttributesBuilder = {
    val runtimeConfig = batchConfiguration.runtimeConfig
    StandardValidatedRuntimeAttributesBuilder.default(runtimeConfig).withValidation(
      cpuValidation(runtimeConfig),
      cpuPlatformValidation(runtimeConfig),
      noAddressValidation(runtimeConfig),
      zonesValidation(runtimeConfig),
      dockerValidation
    )
  }

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes, runtimeAttrsConfig: Option[Config]): GcpBatchRuntimeAttributes = {
    val cpu: Int Refined Positive = RuntimeAttributesValidation.extract(cpuValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val cpuPlatform: Option[String] = RuntimeAttributesValidation.extractOption(cpuPlatformValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)
    val docker: String = RuntimeAttributesValidation.extract(dockerValidation, validatedRuntimeAttributes)
    val failOnStderr: Boolean = RuntimeAttributesValidation.extract(failOnStderrValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val noAddress: Boolean = RuntimeAttributesValidation.extract(noAddressValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val zones: Vector[String] = RuntimeAttributesValidation.extract(ZonesValidation, validatedRuntimeAttributes)

    new GcpBatchRuntimeAttributes(
      cpu,
      cpuPlatform,
      zones,
      docker,
      failOnStderr,
      noAddress)
  }


}
