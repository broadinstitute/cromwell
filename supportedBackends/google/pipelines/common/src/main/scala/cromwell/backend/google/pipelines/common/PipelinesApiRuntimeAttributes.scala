package cromwell.backend.google.pipelines.common

import cats.data.Validated._
import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import com.typesafe.config.Config
import common.validation.ErrorOr._
import cromwell.backend.google.pipelines.common.GpuResource.GpuType
import cromwell.backend.google.pipelines.common.io.{PipelinesApiAttachedDisk, PipelinesApiWorkingDisk}
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation.{BooleanRuntimeAttributesValidation, _}
import eu.timepit.refined._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import wdl4s.parser.MemoryUnit
import wom.RuntimeAttributesKeys
import wom.format.MemorySize
import wom.types._
import wom.values._

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

final case class GpuResource(gpuType: GpuType, gpuCount: Int Refined Positive, nvidiaDriverVersion: String = GpuResource.DefaultNvidiaDriverVersion)

final case class PipelinesApiRuntimeAttributes(cpu: Int Refined Positive,
                                               cpuPlatform: Option[String],
                                               gpuResource: Option[GpuResource],
                                               zones: Vector[String],
                                               preemptible: Int,
                                               bootDiskSize: Int,
                                               memory: MemorySize,
                                               disks: Seq[PipelinesApiAttachedDisk],
                                               dockerImage: String,
                                               failOnStderr: Boolean,
                                               continueOnReturnCode: ContinueOnReturnCode,
                                               noAddress: Boolean)

object PipelinesApiRuntimeAttributes {

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
  private val DisksDefaultValue = WomString(s"${PipelinesApiWorkingDisk.Name} 10 SSD")

  val CpuPlatformKey = "cpuPlatform"
  private val cpuPlatformValidationInstance = new StringRuntimeAttributesValidation(CpuPlatformKey).optional

  private val MemoryDefaultValue = "2048 MB"

  private def cpuValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int Refined Positive] = CpuValidation.instance
    .withDefault(CpuValidation.configDefaultWomValue(runtimeConfig) getOrElse CpuValidation.defaultMin)

  private def cpuPlatformValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] = cpuPlatformValidationInstance

  private def cpuMinValidation(runtimeConfig: Option[Config]):RuntimeAttributesValidation[Int Refined Positive] = CpuValidation.instanceMin
    .withDefault(CpuValidation.configDefaultWomValue(runtimeConfig) getOrElse CpuValidation.defaultMin)

  private def gpuTypeValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[GpuType] = GpuTypeValidation.optional

  val GpuDriverVersionKey = "nvidiaDriverVersion"
  private def gpuDriverValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] = new StringRuntimeAttributesValidation(GpuDriverVersionKey).optional

  private def gpuCountValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Int Refined Positive] = GpuValidation.optional

  private def gpuMinValidation(runtimeConfig: Option[Config]):OptionalRuntimeAttributesValidation[Int Refined Positive] = GpuValidation.optionalMin

  private def failOnStderrValidation(runtimeConfig: Option[Config]) = FailOnStderrValidation.default(runtimeConfig)

  private def continueOnReturnCodeValidation(runtimeConfig: Option[Config]) = ContinueOnReturnCodeValidation.default(runtimeConfig)

  private def disksValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Seq[PipelinesApiAttachedDisk]] = DisksValidation
    .withDefault(DisksValidation.configDefaultWomValue(runtimeConfig) getOrElse DisksDefaultValue)

  private def zonesValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Vector[String]] = ZonesValidation
    .withDefault(ZonesValidation.configDefaultWomValue(runtimeConfig) getOrElse ZonesDefaultValue)

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
      MemoryValidation.configDefaultString(RuntimeAttributesKeys.MemoryMinKey, runtimeConfig) getOrElse MemoryDefaultValue)
  }

  private def bootDiskSizeValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int] = bootDiskValidationInstance
    .withDefault(bootDiskValidationInstance.configDefaultWomValue(runtimeConfig) getOrElse BootDiskDefaultValue)

  private def noAddressValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Boolean] = noAddressValidationInstance
    .withDefault(noAddressValidationInstance.configDefaultWomValue(runtimeConfig) getOrElse NoAddressDefaultValue)

  private val dockerValidation: RuntimeAttributesValidation[String] = DockerValidation.instance

  private val outDirMinValidation: OptionalRuntimeAttributesValidation[MemorySize] = {
    InformationValidation.optional(RuntimeAttributesKeys.OutDirMinKey, MemoryUnit.MB, allowZero = true)
  }

  private val tmpDirMinValidation: OptionalRuntimeAttributesValidation[MemorySize] = {
    InformationValidation.optional(RuntimeAttributesKeys.TmpDirMinKey, MemoryUnit.MB, allowZero = true)
  }

  private val inputDirMinValidation: OptionalRuntimeAttributesValidation[MemorySize] = {
    InformationValidation.optional(RuntimeAttributesKeys.DnaNexusInputDirMinKey, MemoryUnit.MB, allowZero = true)
  }

  def runtimeAttributesBuilder(jesConfiguration: PipelinesApiConfiguration): StandardValidatedRuntimeAttributesBuilder = {
    val runtimeConfig = jesConfiguration.runtimeConfig
    StandardValidatedRuntimeAttributesBuilder.default(runtimeConfig).withValidation(
      gpuCountValidation(runtimeConfig),
      gpuTypeValidation(runtimeConfig),
      gpuDriverValidation(runtimeConfig),
      cpuValidation(runtimeConfig),
      cpuMinValidation(runtimeConfig),
      gpuMinValidation(runtimeConfig),
      disksValidation(runtimeConfig),
      zonesValidation(runtimeConfig),
      preemptibleValidation(runtimeConfig),
      memoryValidation(runtimeConfig),
      memoryMinValidation(runtimeConfig),
      bootDiskSizeValidation(runtimeConfig),
      noAddressValidation(runtimeConfig),
      cpuPlatformValidation(runtimeConfig),
      dockerValidation,
      outDirMinValidation,
      tmpDirMinValidation,
      inputDirMinValidation
    )
  }

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes, runtimeAttrsConfig: Option[Config]): PipelinesApiRuntimeAttributes = {
    val cpu: Int Refined Positive = RuntimeAttributesValidation.extract(cpuValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val cpuPlatform: Option[String] = RuntimeAttributesValidation.extractOption(cpuPlatformValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)

    // GPU
    lazy val gpuType: Option[GpuType] = RuntimeAttributesValidation.extractOption(gpuTypeValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)
    lazy val gpuCount: Option[Int Refined Positive] = RuntimeAttributesValidation.extractOption(gpuCountValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)
    lazy val gpuDriver: Option[String] = RuntimeAttributesValidation.extractOption(gpuDriverValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)

    val gpuResource: Option[GpuResource] = if (gpuType.isDefined || gpuCount.isDefined || gpuDriver.isDefined) {
      Option(GpuResource(gpuType.getOrElse(GpuType.DefaultGpuType), gpuCount.getOrElse(GpuType.DefaultGpuCount), gpuDriver.getOrElse(GpuResource.DefaultNvidiaDriverVersion)))
    } else {
      None
    }

    val zones: Vector[String] = RuntimeAttributesValidation.extract(ZonesValidation, validatedRuntimeAttributes)
    val preemptible: Int = RuntimeAttributesValidation.extract(preemptibleValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val bootDiskSize: Int = RuntimeAttributesValidation.extract(bootDiskSizeValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val memory: MemorySize = RuntimeAttributesValidation.extract(memoryValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val disks: Seq[PipelinesApiAttachedDisk] = RuntimeAttributesValidation.extract(disksValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val docker: String = RuntimeAttributesValidation.extract(dockerValidation, validatedRuntimeAttributes)
    val failOnStderr: Boolean = RuntimeAttributesValidation.extract(failOnStderrValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val continueOnReturnCode: ContinueOnReturnCode = RuntimeAttributesValidation.extract(continueOnReturnCodeValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val noAddress: Boolean = RuntimeAttributesValidation.extract(noAddressValidation(runtimeAttrsConfig), validatedRuntimeAttributes)

    val outDirMin: Option[MemorySize] = RuntimeAttributesValidation.extractOption(outDirMinValidation.key, validatedRuntimeAttributes)
    val tmpDirMin: Option[MemorySize] = RuntimeAttributesValidation.extractOption(tmpDirMinValidation.key, validatedRuntimeAttributes)
    val inputDirMin: Option[MemorySize] = RuntimeAttributesValidation.extractOption(inputDirMinValidation.key, validatedRuntimeAttributes)

    val totalExecutionDiskSizeBytes = List(inputDirMin.map(_.bytes), outDirMin.map(_.bytes), tmpDirMin.map(_.bytes)).flatten.fold(MemorySize(0, MemoryUnit.Bytes).bytes)(_ + _)
    val totalExecutionDiskSize = MemorySize(totalExecutionDiskSizeBytes, MemoryUnit.Bytes)

    val adjustedDisks = disks.adjustWorkingDiskWithNewMin(totalExecutionDiskSize, ())

    new PipelinesApiRuntimeAttributes(
      cpu,
      cpuPlatform,
      gpuResource,
      zones,
      preemptible,
      bootDiskSize,
      memory,
      adjustedDisks,
      docker,
      failOnStderr,
      continueOnReturnCode,
      noAddress
    )
  }
}

object ZonesValidation extends RuntimeAttributesValidation[Vector[String]] {
  override def key: String = PipelinesApiRuntimeAttributes.ZonesKey

  override def coercion: Traversable[WomType] = Set(WomStringType, WomArrayType(WomStringType))

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Vector[String]]] = {
    case WomString(s) => s.split("\\s+").toVector.validNel
    case WomArray(womType, value) if womType.memberType == WomStringType =>
      value.map(_.valueString).toVector.validNel
  }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be either a whitespace separated String or an Array[String]"
}

object DisksValidation extends RuntimeAttributesValidation[Seq[PipelinesApiAttachedDisk]] {
  override def key: String = PipelinesApiRuntimeAttributes.DisksKey

  override def coercion: Traversable[WomType] = Set(WomStringType, WomArrayType(WomStringType))

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Seq[PipelinesApiAttachedDisk]]] = {
    case WomString(value) => validateLocalDisks(value.split(",\\s*").toSeq)
    case WomArray(womType, values) if womType.memberType == WomStringType =>
      validateLocalDisks(values.map(_.valueString))
  }

  private def validateLocalDisks(disks: Seq[String]): ErrorOr[Seq[PipelinesApiAttachedDisk]] = {
    val diskNels: ErrorOr[Seq[PipelinesApiAttachedDisk]] = disks.toList.traverse[ErrorOr, PipelinesApiAttachedDisk](validateLocalDisk)
    val defaulted: ErrorOr[Seq[PipelinesApiAttachedDisk]] = addDefault(diskNels)
    defaulted
  }

  private def validateLocalDisk(disk: String): ErrorOr[PipelinesApiAttachedDisk] = {
    PipelinesApiAttachedDisk.parse(disk) match {
      case scala.util.Success(attachedDisk) => attachedDisk.validNel
      case scala.util.Failure(ex) => ex.getMessage.invalidNel
    }
  }

  private def addDefault(disksNel: ErrorOr[Seq[PipelinesApiAttachedDisk]]): ErrorOr[Seq[PipelinesApiAttachedDisk]] = {
    disksNel map {
      case disks if disks.exists(_.name == PipelinesApiWorkingDisk.Name) => disks
      case disks => disks :+ PipelinesApiWorkingDisk.Default
    }
  }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be a comma separated String or Array[String]"
}
