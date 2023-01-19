package cromwell.backend.google.pipelines.batch

import com.typesafe.config.Config
//import cromwell.backend.google.pipelines.common.PipelinesApiRuntimeAttributes.dockerValidation
import cromwell.backend.validation.{DockerValidation, StringRuntimeAttributesValidation, ValidatedRuntimeAttributes}
//import cromwell.backend.google.pipelines.common.PipelinesApiConfiguration
//import cromwell.backend.google.pipelines.common.PipelinesApiRuntimeAttributes.cpuValidation
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation.{OptionalRuntimeAttributesValidation, RuntimeAttributesValidation}


case class GcpBatchRuntimeAttributes(cpuPlatform: Option[String], dockerImage: String)

object GcpBatchRuntimeAttributes {


  val CpuPlatformKey = "cpuPlatform"
  private val cpuPlatformValidationInstance = new StringRuntimeAttributesValidation(CpuPlatformKey)
    .optional
  // via `gcloud compute zones describe us-central1-a`
  //val CpuPlatformIntelCascadeLakeValue = "Intel Cascade Lake"
  //val CpuPlatformAMDRomeValue = "AMD Rome"

  private def cpuPlatformValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] = cpuPlatformValidationInstance
  private val dockerValidation: RuntimeAttributesValidation[String] = DockerValidation.instance

  def runtimeAttributesBuilder(batchConfiguration: GcpBatchConfiguration): StandardValidatedRuntimeAttributesBuilder = {
    val runtimeConfig = batchConfiguration.runtimeConfig
    StandardValidatedRuntimeAttributesBuilder.default(runtimeConfig).withValidation(
      cpuPlatformValidation(runtimeConfig),
      dockerValidation
    )
  }

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes, runtimeAttrsConfig: Option[Config]): GcpBatchRuntimeAttributes = {
    val cpuPlatform: Option[String] = RuntimeAttributesValidation.extractOption(cpuPlatformValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)
    val docker: String = RuntimeAttributesValidation.extract(dockerValidation, validatedRuntimeAttributes)

    new GcpBatchRuntimeAttributes(cpuPlatform, docker)
  }


}
