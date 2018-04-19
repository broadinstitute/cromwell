package cromwell.backend.google.pipelines.common

import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
import com.typesafe.config.Config
import common.validation.ErrorOr._
import cromwell.backend.google.pipelines.common.GpuResource.GpuType.GpuType
import cromwell.backend.google.pipelines.common.io.{JesAttachedDisk, JesWorkingDisk}
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation.{BooleanRuntimeAttributesValidation, _}
import wom.RuntimeAttributesKeys
import wom.format.MemorySize
import wom.types._
import wom.values._

object GpuResource {
  object GpuType extends Enumeration {
    type GpuType = Value
    val NVIDIATeslaP100 = Value("nvidia-tesla-p100")
    val NVIDIATeslaK80 = Value("nvidia-tesla-k80")
  }
}

final case class GpuResource(gpuType: GpuType, gpuCount: Int)

final case class PipelinesApiRuntimeAttributes(cpu: Int,
                                               gpuResource: Option[GpuResource],
                                               zones: Vector[String],
                                               preemptible: Int,
                                               bootDiskSize: Int,
                                               memory: MemorySize,
                                               disks: Seq[JesAttachedDisk],
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
  private val DisksDefaultValue = WomString(s"${JesWorkingDisk.Name} 10 SSD")

  private val MemoryDefaultValue = "2 GB"

  private def cpuValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int] = CpuValidation.instance
    .withDefault(CpuValidation.configDefaultWomValue(runtimeConfig) getOrElse CpuValidation.defaultMin)

  private def cpuMinValidation(runtimeConfig: Option[Config]):RuntimeAttributesValidation[Int] = CpuValidation.instanceMin
    .withDefault(CpuValidation.configDefaultWomValue(runtimeConfig) getOrElse CpuValidation.defaultMin)

  private def gpuTypeValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[GpuType] = GpuTypeValidation.optional

  private def gpuValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Int] = GpuValidation.optional

  private def gpuMinValidation(runtimeConfig: Option[Config]):OptionalRuntimeAttributesValidation[Int] = GpuValidation.optionalMin

  private def failOnStderrValidation(runtimeConfig: Option[Config]) = FailOnStderrValidation.default(runtimeConfig)

  private def continueOnReturnCodeValidation(runtimeConfig: Option[Config]) = ContinueOnReturnCodeValidation.default(runtimeConfig)

  private def disksValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Seq[JesAttachedDisk]] = DisksValidation
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

  def runtimeAttributesBuilder(jesConfiguration: PipelinesApiConfiguration): StandardValidatedRuntimeAttributesBuilder = {
    val runtimeConfig = jesConfiguration.runtimeConfig
    StandardValidatedRuntimeAttributesBuilder.default(runtimeConfig).withValidation(
      gpuValidation(runtimeConfig),
      gpuTypeValidation(runtimeConfig),
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
      dockerValidation
    )
  }

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes, runtimeAttrsConfig: Option[Config]): PipelinesApiRuntimeAttributes = {
    val cpu: Int = RuntimeAttributesValidation.extract(cpuValidation(runtimeAttrsConfig), validatedRuntimeAttributes)

    // GPU
    val gpuType: Option[GpuType] = RuntimeAttributesValidation.extractOption(gpuTypeValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)
    val gpu: Option[Int] = RuntimeAttributesValidation.extractOption(gpuValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)
    val gpuResource = (gpuType, gpu) match {
      case (Some(t), Some(g)) => Option(GpuResource(t, g))
      case (Some(_), None) => throw new RuntimeException(s"Please specify how many GPU should be attached to the instance.")
      case (None, Some(_)) => throw new RuntimeException(s"Please specify a GPU type: ${GpuResource.GpuType.values.mkString(", ")}")
      case (None, None) => None
    }

    val zones: Vector[String] = RuntimeAttributesValidation.extract(ZonesValidation, validatedRuntimeAttributes)
    val preemptible: Int = RuntimeAttributesValidation.extract(preemptibleValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val bootDiskSize: Int = RuntimeAttributesValidation.extract(bootDiskSizeValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val memory: MemorySize = RuntimeAttributesValidation.extract(memoryValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val disks: Seq[JesAttachedDisk] = RuntimeAttributesValidation.extract(disksValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val docker: String = RuntimeAttributesValidation.extract(dockerValidation, validatedRuntimeAttributes)
    val failOnStderr: Boolean = RuntimeAttributesValidation.extract(failOnStderrValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val continueOnReturnCode: ContinueOnReturnCode = RuntimeAttributesValidation.extract(continueOnReturnCodeValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val noAddress: Boolean = RuntimeAttributesValidation.extract(noAddressValidation(runtimeAttrsConfig), validatedRuntimeAttributes)

    new PipelinesApiRuntimeAttributes(
      cpu,
      gpuResource,
      zones,
      preemptible,
      bootDiskSize,
      memory,
      disks,
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

object DisksValidation extends RuntimeAttributesValidation[Seq[JesAttachedDisk]] {
  override def key: String = PipelinesApiRuntimeAttributes.DisksKey

  override def coercion: Traversable[WomType] = Set(WomStringType, WomArrayType(WomStringType))

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Seq[JesAttachedDisk]]] = {
    case WomString(value) => validateLocalDisks(value.split(",\\s*").toSeq)
    case WomArray(womType, values) if womType.memberType == WomStringType =>
      validateLocalDisks(values.map(_.valueString))
  }

  private def validateLocalDisks(disks: Seq[String]): ErrorOr[Seq[JesAttachedDisk]] = {
    val diskNels: Seq[ErrorOr[JesAttachedDisk]] = disks map validateLocalDisk
    val sequenced: ErrorOr[Seq[JesAttachedDisk]] = sequenceNels(diskNels)
    val defaulted: ErrorOr[Seq[JesAttachedDisk]] = addDefault(sequenced)
    defaulted
  }

  private def validateLocalDisk(disk: String): ErrorOr[JesAttachedDisk] = {
    JesAttachedDisk.parse(disk) match {
      case scala.util.Success(attachedDisk) => attachedDisk.validNel
      case scala.util.Failure(ex) => ex.getMessage.invalidNel
    }
  }

  private def sequenceNels(nels: Seq[ErrorOr[JesAttachedDisk]]): ErrorOr[Seq[JesAttachedDisk]] = {
    val emptyDiskNel: ErrorOr[Vector[JesAttachedDisk]] = Vector.empty[JesAttachedDisk].validNel
    val disksNel: ErrorOr[Vector[JesAttachedDisk]] = nels.foldLeft(emptyDiskNel) {
      (acc, v) => (acc, v) mapN { (a, v) => a :+ v }
    }
    disksNel
  }

  private def addDefault(disksNel: ErrorOr[Seq[JesAttachedDisk]]): ErrorOr[Seq[JesAttachedDisk]] = {
    disksNel map {
      case disks if disks.exists(_.name == JesWorkingDisk.Name) => disks
      case disks => disks :+ JesWorkingDisk.Default
    }
  }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be a comma separated String or Array[String]"
}
