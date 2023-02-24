package cromwell.backend.google.pipelines.batch

import com.typesafe.config.Config
//import cromwell.backend.google.pipelines.batch.GpuResource.GpuType
import cromwell.backend.google.pipelines.common.ZonesValidation
//import cromwell.backend.google.pipelines.common.{GpuTypeValidation, GpuValidation, ZonesValidation}
//import cromwell.backend.google.pipelines.common.{GpuTypeValidation, GpuValidation, ZonesValidation}
import cromwell.backend.validation.BooleanRuntimeAttributesValidation
import wom.values.{WomBoolean, WomInteger, WomString}
import cromwell.backend.validation.{DockerValidation, StringRuntimeAttributesValidation, ValidatedRuntimeAttributes}
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation.{OptionalRuntimeAttributesValidation, RuntimeAttributesValidation}
import eu.timepit.refined.api.Refined
import cromwell.backend.validation._
import eu.timepit.refined.numeric.Positive
//import eu.timepit.refined.refineMV
import wom.RuntimeAttributesKeys
import wom.format.MemorySize
/*
object GpuResource {
  val DefaultNvidiaDriverVersion = "418.87.00"
  final case class GpuType(name: String) {
    override def toString: String = name
  }
  object GpuType {
    val NVIDIATeslaP100 = GpuType("nvidia-tesla-p100")
    val NVIDIATeslaK80 = GpuType("nvidia-tesla-k80")

    val DefaultGpuType: GpuType = NVIDIATeslaK80
    val DefaultGpuCount: Int Refined Positive = refineMV[Positive](1)
    val MoreDetailsURL = "https://cloud.google.com/compute/docs/gpus/"
  }
}*/

//final case class GpuResource(gpuType: GpuType, gpuCount: Int Refined Positive, nvidiaDriverVersion: String = GpuResource.DefaultNvidiaDriverVersion)

final case class GcpBatchRuntimeAttributes(
                                      cpu: Int Refined Positive,
                                      cpuPlatform: Option[String],
                                      //gpuResource: Option[GpuResource],
                                      zones: Vector[String],
                                      preemptible: Int,
                                      bootDiskSize: Int,
                                      memory: MemorySize,
                                      dockerImage: String,
                                      failOnStderr: Boolean,
                                      noAddress: Boolean,
                                      useDockerImageCache: Option[Boolean],
                                      checkpointFilename: Option[String])

object GcpBatchRuntimeAttributes {

  val ZonesKey = "zones"
  private val ZonesDefaultValue = WomString("us-central1-b")

  val PreemptibleKey = "preemptible"
  private val preemptibleValidationInstance = new IntRuntimeAttributesValidation(PreemptibleKey)
  private val PreemptibleDefaultValue = WomInteger(0)

  val BootDiskSizeKey = "bootDiskSizeGb"
  private val bootDiskValidationInstance = new IntRuntimeAttributesValidation(BootDiskSizeKey)
  private val BootDiskDefaultValue = WomInteger(10)

  val NoAddressKey = "noAddress"
  private val noAddressValidationInstance = new BooleanRuntimeAttributesValidation(NoAddressKey)
  private val NoAddressDefaultValue = WomBoolean(false)


  val CpuPlatformKey = "cpuPlatform"
  private val cpuPlatformValidationInstance = new StringRuntimeAttributesValidation(CpuPlatformKey)
    .optional
  // via `gcloud compute zones describe us-central1-a`
  val CpuPlatformIntelCascadeLakeValue = "Intel Cascade Lake"
  val CpuPlatformAMDRomeValue = "AMD Rome"

  private def cpuMinValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int Refined Positive] = CpuValidation
    .instanceMin
    .withDefault(CpuValidation.configDefaultWomValue(runtimeConfig) getOrElse CpuValidation.defaultMin)

  val UseDockerImageCacheKey = "useDockerImageCache"
  private val useDockerImageCacheValidationInstance = new BooleanRuntimeAttributesValidation(UseDockerImageCacheKey)
    .optional

  val CheckpointFileKey = "checkpointFile"
  private val checkpointFileValidationInstance = new StringRuntimeAttributesValidation(CheckpointFileKey).optional

  private val MemoryDefaultValue = "2048 MB"

  private def cpuValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int Refined Positive] = CpuValidation
    .instance
    .withDefault(CpuValidation
      .configDefaultWomValue(runtimeConfig) getOrElse CpuValidation
      .defaultMin)
  private def cpuPlatformValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] = cpuPlatformValidationInstance
  //private def gpuTypeValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[GpuType] = GpuTypeValidation.optional

  val GpuDriverVersionKey = "nvidiaDriverVersion"
  //private def gpuDriverValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] = new StringRuntimeAttributesValidation(GpuDriverVersionKey).optional
  //private def gpuCountValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Int Refined Positive] = GpuValidation.optional
  //private def gpuMinValidation(runtimeConfig: Option[Config]):OptionalRuntimeAttributesValidation[Int Refined Positive] = GpuValidation.optionalMin

  private val dockerValidation: RuntimeAttributesValidation[String] = DockerValidation.instance



  private def failOnStderrValidation(runtimeConfig: Option[Config]) = FailOnStderrValidation.default(runtimeConfig)
  private def zonesValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Vector[String]] = ZonesValidation
    .withDefault(ZonesValidation
      .configDefaultWomValue(runtimeConfig) getOrElse ZonesDefaultValue)

  private def preemptibleValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int] = preemptibleValidationInstance
    .withDefault(preemptibleValidationInstance.configDefaultWomValue(runtimeConfig) getOrElse PreemptibleDefaultValue)

  private def memoryValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[MemorySize] = {
    MemoryValidation.withDefaultMemory(
      RuntimeAttributesKeys.MemoryKey,
      MemoryValidation.configDefaultString(RuntimeAttributesKeys.MemoryKey, runtimeConfig) getOrElse MemoryDefaultValue)
  }

  private def memoryMinValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[MemorySize] = {
    MemoryValidation.withDefaultMemory(
      RuntimeAttributesKeys.MemoryMinKey,
      MemoryValidation
        .configDefaultString(RuntimeAttributesKeys.MemoryMinKey, runtimeConfig) getOrElse MemoryDefaultValue)
  }


  private def bootDiskSizeValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int] = bootDiskValidationInstance
    .withDefault(bootDiskValidationInstance.configDefaultWomValue(runtimeConfig) getOrElse BootDiskDefaultValue)

  private def noAddressValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Boolean] = noAddressValidationInstance
    .withDefault(noAddressValidationInstance
      .configDefaultWomValue(runtimeConfig) getOrElse NoAddressDefaultValue)

  private def useDockerImageCacheValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Boolean] =
    useDockerImageCacheValidationInstance

  def runtimeAttributesBuilder(batchConfiguration: GcpBatchConfiguration): StandardValidatedRuntimeAttributesBuilder = {
    val runtimeConfig = batchConfiguration.runtimeConfig
    StandardValidatedRuntimeAttributesBuilder.default(runtimeConfig).withValidation(
      //gpuCountValidation(runtimeConfig),
      //gpuTypeValidation(runtimeConfig),
      //gpuDriverValidation(runtimeConfig),
      cpuValidation(runtimeConfig),
      cpuPlatformValidation(runtimeConfig),
      cpuMinValidation(runtimeConfig),
      //gpuMinValidation(runtimeConfig),
      noAddressValidation(runtimeConfig),
      zonesValidation(runtimeConfig),
      preemptibleValidation(runtimeConfig),
      memoryValidation(runtimeConfig),
      memoryMinValidation(runtimeConfig),
      bootDiskSizeValidation(runtimeConfig),
      useDockerImageCacheValidation(runtimeConfig),
      checkpointFileValidationInstance,
      dockerValidation
    )
  }

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes, runtimeAttrsConfig: Option[Config]): GcpBatchRuntimeAttributes = {
    val cpu: Int Refined Positive = RuntimeAttributesValidation.extract(cpuValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val cpuPlatform: Option[String] = RuntimeAttributesValidation.extractOption(cpuPlatformValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)
    val checkpointFileName: Option[String] = RuntimeAttributesValidation.extractOption(checkpointFileValidationInstance.key, validatedRuntimeAttributes)

    // GPU
    //lazy val gpuType: Option[GpuType] = RuntimeAttributesValidation
    //  .extractOption(gpuTypeValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)
    //lazy val gpuCount: Option[Int Refined Positive] = RuntimeAttributesValidation
      //.extractOption(gpuCountValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)
    //lazy val gpuDriver: Option[String] = RuntimeAttributesValidation
      //.extractOption(gpuDriverValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)

    //val gpuResource: Option[GpuResource] = if (gpuType.isDefined || gpuCount.isDefined || gpuDriver.isDefined) {
    //  Option(GpuResource(gpuType.getOrElse(GpuType.DefaultGpuType), gpuCount
    //    .getOrElse(GpuType.DefaultGpuCount), gpuDriver.getOrElse(GpuResource.DefaultNvidiaDriverVersion)))
    //} else {
    //  None
    //}

    val docker: String = RuntimeAttributesValidation.extract(dockerValidation, validatedRuntimeAttributes)
    val failOnStderr: Boolean = RuntimeAttributesValidation.extract(failOnStderrValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val noAddress: Boolean = RuntimeAttributesValidation.extract(noAddressValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val zones: Vector[String] = RuntimeAttributesValidation.extract(ZonesValidation, validatedRuntimeAttributes)
    val preemptible: Int = RuntimeAttributesValidation.extract(preemptibleValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val bootDiskSize: Int = RuntimeAttributesValidation.extract(bootDiskSizeValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val memory: MemorySize = RuntimeAttributesValidation.extract(memoryValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val useDockerImageCache: Option[Boolean] = RuntimeAttributesValidation.extractOption(useDockerImageCacheValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)

    new GcpBatchRuntimeAttributes(
      cpu,
      cpuPlatform,
      //gpuResource,
      zones,
      preemptible,
      bootDiskSize,
      memory,
      docker,
      failOnStderr,
      noAddress,
      useDockerImageCache,
      checkpointFileName)
  }


}
