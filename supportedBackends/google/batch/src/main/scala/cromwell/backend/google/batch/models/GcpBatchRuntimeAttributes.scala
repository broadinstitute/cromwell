package cromwell.backend.google.batch.models

import cats.implicits.{catsSyntaxValidatedId, toTraverseOps}
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.google.batch.io.{GcpBatchAttachedDisk, GcpBatchWorkingDisk}
import cromwell.backend.google.batch.models.GpuResource.GpuType
import cromwell.backend.google.batch.util.{GpuTypeValidation, GpuValidation}
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation._
import eu.timepit.refined._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import wom.RuntimeAttributesKeys
import wom.format.MemorySize
import wom.types.{WomArrayType, WomStringType, WomType}
import wom.values.{WomArray, WomBoolean, WomInteger, WomString, WomValue}

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
}

final case class GpuResource(gpuType: GpuType, gpuCount: Int Refined Positive)

final case class GcpBatchRuntimeAttributes(cpu: Int Refined Positive,
                                           cpuPlatform: Option[String],
                                           gpuResource: Option[GpuResource],
                                           zones: Vector[String],
                                           preemptible: Int,
                                           bootDiskSize: Int,
                                           memory: MemorySize,
                                           disks: Seq[GcpBatchAttachedDisk],
                                           dockerImage: String,
                                           failOnStderr: Boolean,
                                           continueOnReturnCode: ContinueOnReturnCode,
                                           noAddress: Boolean,
                                           useDockerImageCache: Option[Boolean],
                                           checkpointFilename: Option[String]
)

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

  val DisksKey = "disks"
  private val DisksDefaultValue = WomString(s"${GcpBatchWorkingDisk.Name} 10 SSD")

  val CpuPlatformKey = "cpuPlatform"
  private val cpuPlatformValidationInstance = new StringRuntimeAttributesValidation(CpuPlatformKey).optional
  // via `gcloud compute zones describe us-central1-a`
  val CpuPlatformIntelCascadeLakeValue = "Intel Cascade Lake"
  val CpuPlatformAMDRomeValue = "AMD Rome"

  val UseDockerImageCacheKey = "useDockerImageCache"
  private val useDockerImageCacheValidationInstance = new BooleanRuntimeAttributesValidation(
    UseDockerImageCacheKey
  ).optional

  val CheckpointFileKey = "checkpointFile"
  private val checkpointFileValidationInstance = new StringRuntimeAttributesValidation(CheckpointFileKey).optional

  private val MemoryDefaultValue = "2048 MB"

  private def cpuValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int Refined Positive] =
    CpuValidation.instance
      .withDefault(
        CpuValidation
          .configDefaultWomValue(runtimeConfig) getOrElse CpuValidation.defaultMin
      )
  private def cpuPlatformValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] =
    cpuPlatformValidationInstance
  private def gpuTypeValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[GpuType] =
    GpuTypeValidation.optional

  val GpuDriverVersionKey = "nvidiaDriverVersion"
  private def gpuDriverValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] =
    new StringRuntimeAttributesValidation(GpuDriverVersionKey).optional

  private def gpuCountValidation(
    runtimeConfig: Option[Config]
  ): OptionalRuntimeAttributesValidation[Int Refined Positive] = GpuValidation.optional

  private val dockerValidation: RuntimeAttributesValidation[String] = DockerValidation.instance

  private def failOnStderrValidation(runtimeConfig: Option[Config]) = FailOnStderrValidation.default(runtimeConfig)

  private def continueOnReturnCodeValidation(runtimeConfig: Option[Config]) =
    ContinueOnReturnCodeValidation.default(runtimeConfig)

  private def disksValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Seq[GcpBatchAttachedDisk]] =
    DisksValidation
      .withDefault(DisksValidation.configDefaultWomValue(runtimeConfig) getOrElse DisksDefaultValue)

  private def zonesValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Vector[String]] =
    ZonesValidation
      .withDefault(
        ZonesValidation
          .configDefaultWomValue(runtimeConfig) getOrElse ZonesDefaultValue
      )

  private def preemptibleValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int] =
    preemptibleValidationInstance
      .withDefault(preemptibleValidationInstance.configDefaultWomValue(runtimeConfig) getOrElse PreemptibleDefaultValue)

  private def memoryValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[MemorySize] =
    MemoryValidation.withDefaultMemory(
      RuntimeAttributesKeys.MemoryKey,
      MemoryValidation.configDefaultString(RuntimeAttributesKeys.MemoryKey, runtimeConfig) getOrElse MemoryDefaultValue
    )

  private def bootDiskSizeValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int] =
    bootDiskValidationInstance
      .withDefault(bootDiskValidationInstance.configDefaultWomValue(runtimeConfig) getOrElse BootDiskDefaultValue)

  private def noAddressValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Boolean] =
    noAddressValidationInstance
      .withDefault(
        noAddressValidationInstance
          .configDefaultWomValue(runtimeConfig) getOrElse NoAddressDefaultValue
      )

  private def useDockerImageCacheValidation(
    runtimeConfig: Option[Config]
  ): OptionalRuntimeAttributesValidation[Boolean] =
    useDockerImageCacheValidationInstance

  def runtimeAttributesBuilder(batchConfiguration: GcpBatchConfiguration): StandardValidatedRuntimeAttributesBuilder = {
    val runtimeConfig = batchConfiguration.runtimeConfig
    StandardValidatedRuntimeAttributesBuilder
      .default(runtimeConfig)
      .withValidation(
        gpuCountValidation(runtimeConfig),
        gpuTypeValidation(runtimeConfig),
        gpuDriverValidation(runtimeConfig),
        cpuValidation(runtimeConfig),
        cpuPlatformValidation(runtimeConfig),
        disksValidation(runtimeConfig),
        noAddressValidation(runtimeConfig),
        zonesValidation(runtimeConfig),
        preemptibleValidation(runtimeConfig),
        memoryValidation(runtimeConfig),
        bootDiskSizeValidation(runtimeConfig),
        useDockerImageCacheValidation(runtimeConfig),
        checkpointFileValidationInstance,
        dockerValidation
      )
  }

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes,
            runtimeAttrsConfig: Option[Config]
  ): GcpBatchRuntimeAttributes = {
    val cpu: Int Refined Positive =
      RuntimeAttributesValidation.extract(cpuValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val cpuPlatform: Option[String] = RuntimeAttributesValidation.extractOption(
      cpuPlatformValidation(runtimeAttrsConfig).key,
      validatedRuntimeAttributes
    )
    val checkpointFileName: Option[String] =
      RuntimeAttributesValidation.extractOption(checkpointFileValidationInstance.key, validatedRuntimeAttributes)

    // GPU
    lazy val gpuType: Option[GpuType] = RuntimeAttributesValidation
      .extractOption(gpuTypeValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)
    lazy val gpuCount: Option[Int Refined Positive] = RuntimeAttributesValidation
      .extractOption(gpuCountValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)
    lazy val gpuDriver: Option[String] =
      RuntimeAttributesValidation.extractOption(gpuDriverValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)

    val gpuResource: Option[GpuResource] = if (gpuType.isDefined || gpuCount.isDefined || gpuDriver.isDefined) {
      Option(
        GpuResource(gpuType.getOrElse(GpuType.DefaultGpuType),
                    gpuCount
                      .getOrElse(GpuType.DefaultGpuCount)
        )
      )
    } else {
      None
    }

    val docker: String = RuntimeAttributesValidation.extract(dockerValidation, validatedRuntimeAttributes)
    val failOnStderr: Boolean =
      RuntimeAttributesValidation.extract(failOnStderrValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val continueOnReturnCode: ContinueOnReturnCode = RuntimeAttributesValidation.extract(
      continueOnReturnCodeValidation(runtimeAttrsConfig),
      validatedRuntimeAttributes
    )
    val noAddress: Boolean =
      RuntimeAttributesValidation.extract(noAddressValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val zones: Vector[String] = RuntimeAttributesValidation.extract(ZonesValidation, validatedRuntimeAttributes)
    val preemptible: Int =
      RuntimeAttributesValidation.extract(preemptibleValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val bootDiskSize: Int =
      RuntimeAttributesValidation.extract(bootDiskSizeValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val memory: MemorySize =
      RuntimeAttributesValidation.extract(memoryValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val disks: Seq[GcpBatchAttachedDisk] =
      RuntimeAttributesValidation.extract(disksValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val useDockerImageCache: Option[Boolean] = RuntimeAttributesValidation.extractOption(
      useDockerImageCacheValidation(runtimeAttrsConfig).key,
      validatedRuntimeAttributes
    )

    new GcpBatchRuntimeAttributes(
      cpu,
      cpuPlatform,
      gpuResource,
      zones,
      preemptible,
      bootDiskSize,
      memory,
      disks,
      docker,
      failOnStderr,
      continueOnReturnCode,
      noAddress,
      useDockerImageCache,
      checkpointFileName
    )
  }

}

object ZonesValidation extends RuntimeAttributesValidation[Vector[String]] {
  override def key: String = GcpBatchRuntimeAttributes.ZonesKey

  override def coercion: Iterable[WomType] = Set(WomStringType, WomArrayType(WomStringType))

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Vector[String]]] = {
    case WomString(s) => s.split("\\s+").toVector.validNel
    case WomArray(womType, value) if womType.memberType == WomStringType =>
      value.map(_.valueString).toVector.validNel
  }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be either a whitespace separated String or an Array[String]"
}

object DisksValidation extends RuntimeAttributesValidation[Seq[GcpBatchAttachedDisk]] {
  override def key: String = GcpBatchRuntimeAttributes.DisksKey

  override def coercion: Iterable[WomType] = Set(WomStringType, WomArrayType(WomStringType))

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Seq[GcpBatchAttachedDisk]]] = {
    case WomString(value) => validateLocalDisks(value.split(",\\s*").toSeq)
    case WomArray(womType, values) if womType.memberType == WomStringType =>
      validateLocalDisks(values.map(_.valueString))
  }

  private def validateLocalDisks(disks: Seq[String]): ErrorOr[Seq[GcpBatchAttachedDisk]] = {
    val diskNels: ErrorOr[Seq[GcpBatchAttachedDisk]] =
      disks.toList.traverse[ErrorOr, GcpBatchAttachedDisk](validateLocalDisk)
    val defaulted: ErrorOr[Seq[GcpBatchAttachedDisk]] = addDefault(diskNels)
    defaulted
  }

  private def validateLocalDisk(disk: String): ErrorOr[GcpBatchAttachedDisk] =
    GcpBatchAttachedDisk.parse(disk) match {
      case scala.util.Success(attachedDisk) => attachedDisk.validNel
      case scala.util.Failure(ex) => ex.getMessage.invalidNel
    }

  private def addDefault(disksNel: ErrorOr[Seq[GcpBatchAttachedDisk]]): ErrorOr[Seq[GcpBatchAttachedDisk]] =
    disksNel map {
      case disks if disks.exists(_.name == GcpBatchWorkingDisk.Name) => disks
      case disks => disks :+ GcpBatchWorkingDisk.Default
    }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be a comma separated String or Array[String]"
}
