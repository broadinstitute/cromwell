package cromwell.backend.google.batch.models

import cats.implicits.{catsSyntaxValidatedId, toTraverseOps}
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.google.batch.io.{GcpBatchAttachedDisk, GcpBatchWorkingDisk}
import cromwell.backend.google.batch.models.GcpBatchRuntimeAttributes.BootDiskSizeKey
import cromwell.backend.google.batch.models.GpuResource.GpuType
import cromwell.backend.google.batch.util.{GpuTypeValidation, GpuValidation}
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation._
import eu.timepit.refined._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import wom.RuntimeAttributesKeys
import wom.format.MemorySize
import wom.types.{WomArrayType, WomIntegerType, WomStringType, WomType}
import wom.values.{WomArray, WomBoolean, WomInteger, WomOptionalValue, WomString, WomValue}

object GpuResource {

  final case class GpuType(name: String) {
    override def toString: String = name
  }

  object GpuType {
    val NVIDIATeslaT4 = GpuType("nvidia-tesla-t4")

    val DefaultGpuType: GpuType = NVIDIATeslaT4
    val DefaultGpuCount: Int Refined Positive = refineMV[Positive](1)
    val MoreDetailsURL = "https://cloud.google.com/compute/docs/gpus/"
  }
}

final case class GpuResource(gpuType: GpuType, gpuCount: Int Refined Positive)

final case class Machine(machineType: String)

final case class GcpBatchRuntimeAttributes(cpu: Int Refined Positive,
                                           cpuPlatform: Option[String],
                                           gpuResource: Option[GpuResource],
                                           zones: Vector[String],
                                           preemptible: Int,
                                           bootDiskSize: Int,
                                           memory: MemorySize,
                                           machine: Option[Machine] = None,
                                           disks: Seq[GcpBatchAttachedDisk],
                                           dockerImage: String,
                                           failOnStderr: Boolean,
                                           continueOnReturnCode: ContinueOnReturnCode,
                                           noAddress: Boolean,
                                           checkpointFilename: Option[String]
)

object GcpBatchRuntimeAttributes {

  val ZonesKey = "zones"
  private val ZonesDefaultValue = WomString("us-central1-b")

  val PreemptibleKey = "preemptible"
  private val preemptibleValidationInstance = new IntRuntimeAttributesValidation(PreemptibleKey)
  private val PreemptibleDefaultValue = WomInteger(0)

  val BootDiskSizeKey = "bootDiskSizeGb"
  // This is the smallest size we will actually get, see https://cloud.google.com/batch/docs/vm-os-environment-overview#default
  val BootDiskDefaultValue = WomInteger(30)

  val NoAddressKey = "noAddress"
  private val noAddressValidationInstance = new BooleanRuntimeAttributesValidation(NoAddressKey)
  private val NoAddressDefaultValue = WomBoolean(false)

  val DisksKey = "disks"
  private val DisksDefaultValue = WomString(s"${GcpBatchWorkingDisk.Name} 10 SSD")

  val CpuPlatformKey = "cpuPlatform"
  private val cpuPlatformValidationInstance = new StringRuntimeAttributesValidation(CpuPlatformKey).optional
  // via `gcloud compute zones describe us-central1-a`
  val CpuPlatformIntelCascadeLakeValue = "Intel Cascade Lake"
  val CpuPlatformIntelIceLakeValue = "Intel Ice Lake"
  val CpuPlatformAMDRomeValue = "AMD Rome"

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

  private def gpuCountValidation(
    runtimeConfig: Option[Config]
  ): OptionalRuntimeAttributesValidation[Int Refined Positive] = GpuValidation.optional

  // As of WDL 1.1 these two are aliases of each other
  private val dockerValidation: OptionalRuntimeAttributesValidation[Containers] = DockerValidation.instance
  private val containerValidation: OptionalRuntimeAttributesValidation[Containers] = ContainerValidation.instance

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
    new BootDiskSizeValidation(runtimeConfig)

  private def noAddressValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Boolean] =
    noAddressValidationInstance
      .withDefault(
        noAddressValidationInstance
          .configDefaultWomValue(runtimeConfig) getOrElse NoAddressDefaultValue
      )

  def runtimeAttributesBuilder(batchConfiguration: GcpBatchConfiguration): StandardValidatedRuntimeAttributesBuilder = {
    val runtimeConfig = batchConfiguration.runtimeConfig
    StandardValidatedRuntimeAttributesBuilder
      .default(runtimeConfig)
      .withValidation(
        gpuCountValidation(runtimeConfig),
        gpuTypeValidation(runtimeConfig),
        cpuValidation(runtimeConfig),
        cpuPlatformValidation(runtimeConfig),
        disksValidation(runtimeConfig),
        noAddressValidation(runtimeConfig),
        zonesValidation(runtimeConfig),
        preemptibleValidation(runtimeConfig),
        memoryValidation(runtimeConfig),
        bootDiskSizeValidation(runtimeConfig),
        checkpointFileValidationInstance,
        dockerValidation,
        containerValidation
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

    val gpuResource: Option[GpuResource] = if (gpuType.isDefined || gpuCount.isDefined) {
      Option(
        GpuResource(gpuType.getOrElse(GpuType.DefaultGpuType),
                    gpuCount
                      .getOrElse(GpuType.DefaultGpuCount)
        )
      )
    } else {
      None
    }

    val docker: String = Containers.extractContainer(validatedRuntimeAttributes)
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

    new GcpBatchRuntimeAttributes(
      cpu = cpu,
      cpuPlatform = cpuPlatform,
      gpuResource = gpuResource,
      zones = zones,
      preemptible = preemptible,
      bootDiskSize = bootDiskSize,
      memory = memory,
      disks = disks,
      dockerImage = docker,
      failOnStderr = failOnStderr,
      continueOnReturnCode = continueOnReturnCode,
      noAddress = noAddress,
      checkpointFilename = checkpointFileName
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

class BootDiskSizeValidation(runtimeConfig: Option[Config]) extends IntRuntimeAttributesValidation(BootDiskSizeKey) {

  override protected def staticDefaultOption: Option[WomInteger] = Option(WomInteger(0))

  private def defaultBootDiskSize: Int = {
    val configDefault: Option[WomInteger] = this.configDefaultWomValue(runtimeConfig) match {
      case Some(i: WomInteger) => Some(i)
      case _ => None
    }
    configDefault.getOrElse(GcpBatchRuntimeAttributes.BootDiskDefaultValue).value
  }

  // If the user doesn't provide a boot disk size, use the default. If they DO provide a boot disk size,
  // add the default to it to account for the space taken by the user Docker image. See AN-345.
  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Int]] = {
    case WomInteger(value) =>
      (value + defaultBootDiskSize).validNel
    case WomOptionalValue(WomIntegerType, Some(WomInteger(value))) =>
      (value + defaultBootDiskSize).validNel
  }
}
