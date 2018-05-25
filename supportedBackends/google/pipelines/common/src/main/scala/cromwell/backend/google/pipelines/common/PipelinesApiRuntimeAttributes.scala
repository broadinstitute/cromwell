package cromwell.backend.google.pipelines.common

import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
import com.typesafe.config.Config
import common.validation.ErrorOr._
import cromwell.backend.google.pipelines.common.GpuResource.GpuType.GpuType
import cromwell.backend.google.pipelines.common.io.{PipelinesApiAttachedDisk, PipelinesApiWorkingDisk}
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation.{BooleanRuntimeAttributesValidation, _}
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import wom.RuntimeAttributesKeys
import wom.types._
import wom.values._
import net.ceedubs.ficus.Ficus._
import squants.information.Information

object GpuResource {
  val DefaultNvidiaDriverVersion = "390.46"
  object GpuType extends Enumeration {
    type GpuType = Value
    val NVIDIATeslaP100 = Value("nvidia-tesla-p100")
    val NVIDIATeslaK80 = Value("nvidia-tesla-k80")
  }
}

final case class GpuResource(gpuType: GpuType, gpuCount: Int Refined Positive, nvidiaDriverVersion: String = GpuResource.DefaultNvidiaDriverVersion)

final case class PipelinesApiRuntimeAttributes(cpu: Int Refined Positive,
                                               gpuResource: Option[GpuResource],
                                               zones: Vector[String],
                                               preemptible: Int,
                                               bootDiskSize: Int,
                                               memory: Information,
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

  private val MemoryDefaultValue = "2048 MB"

  private def cpuValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int Refined Positive] = CpuValidation.instance
    .withDefault(CpuValidation.configDefaultWomValue(runtimeConfig) getOrElse CpuValidation.defaultMin)

  private def cpuMinValidation(runtimeConfig: Option[Config]):RuntimeAttributesValidation[Int Refined Positive] = CpuValidation.instanceMin
    .withDefault(CpuValidation.configDefaultWomValue(runtimeConfig) getOrElse CpuValidation.defaultMin)

  private def gpuTypeValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[GpuType] = GpuTypeValidation.optional

  private def gpuValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Int Refined Positive] = GpuValidation.optional

  private def gpuMinValidation(runtimeConfig: Option[Config]):OptionalRuntimeAttributesValidation[Int Refined Positive] = GpuValidation.optionalMin

  private def failOnStderrValidation(runtimeConfig: Option[Config]) = FailOnStderrValidation.default(runtimeConfig)

  private def continueOnReturnCodeValidation(runtimeConfig: Option[Config]) = ContinueOnReturnCodeValidation.default(runtimeConfig)

  private def disksValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Seq[PipelinesApiAttachedDisk]] = DisksValidation
    .withDefault(DisksValidation.configDefaultWomValue(runtimeConfig) getOrElse DisksDefaultValue)

  private def zonesValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Vector[String]] = ZonesValidation
    .withDefault(ZonesValidation.configDefaultWomValue(runtimeConfig) getOrElse ZonesDefaultValue)

  private def preemptibleValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int] = preemptibleValidationInstance
    .withDefault(preemptibleValidationInstance.configDefaultWomValue(runtimeConfig) getOrElse PreemptibleDefaultValue)

  private def memoryValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Information] = {
    MemoryValidation.withDefaultMemory(
      RuntimeAttributesKeys.MemoryKey,
      MemoryValidation.configDefaultString(RuntimeAttributesKeys.MemoryKey, runtimeConfig) getOrElse MemoryDefaultValue)
  }

  private def memoryMinValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Information] = {
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
    val cpu: Int Refined Positive = RuntimeAttributesValidation.extract(cpuValidation(runtimeAttrsConfig), validatedRuntimeAttributes)

    // GPU
    lazy val nvidiaDriverVersion = runtimeAttrsConfig.flatMap(_.as[Option[String]]("nvidia-driver-version")).getOrElse(GpuResource.DefaultNvidiaDriverVersion)
    val gpuType: Option[GpuType] = RuntimeAttributesValidation.extractOption(gpuTypeValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)
    val gpu: Option[Int Refined Positive] = RuntimeAttributesValidation.extractOption(gpuValidation(runtimeAttrsConfig).key, validatedRuntimeAttributes)
    val gpuResource = (gpuType, gpu) match {
      case (Some(t), Some(g)) => Option(GpuResource(t, g, nvidiaDriverVersion))
      case (Some(_), None) => throw new RuntimeException(s"Please specify how many GPU should be attached to the instance.")
      case (None, Some(_)) => throw new RuntimeException(s"Please specify a GPU type: ${GpuResource.GpuType.values.mkString(", ")}")
      case (None, None) => None
    }

    val zones: Vector[String] = RuntimeAttributesValidation.extract(ZonesValidation, validatedRuntimeAttributes)
    val preemptible: Int = RuntimeAttributesValidation.extract(preemptibleValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val bootDiskSize: Int = RuntimeAttributesValidation.extract(bootDiskSizeValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val memory: Information = RuntimeAttributesValidation.extract(memoryValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val disks: Seq[PipelinesApiAttachedDisk] = RuntimeAttributesValidation.extract(disksValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
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

object DisksValidation extends RuntimeAttributesValidation[Seq[PipelinesApiAttachedDisk]] {
  override def key: String = PipelinesApiRuntimeAttributes.DisksKey

  override def coercion: Traversable[WomType] = Set(WomStringType, WomArrayType(WomStringType))

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Seq[PipelinesApiAttachedDisk]]] = {
    case WomString(value) => validateLocalDisks(value.split(",\\s*").toSeq)
    case WomArray(womType, values) if womType.memberType == WomStringType =>
      validateLocalDisks(values.map(_.valueString))
  }

  private def validateLocalDisks(disks: Seq[String]): ErrorOr[Seq[PipelinesApiAttachedDisk]] = {
    val diskNels: Seq[ErrorOr[PipelinesApiAttachedDisk]] = disks map validateLocalDisk
    val sequenced: ErrorOr[Seq[PipelinesApiAttachedDisk]] = sequenceNels(diskNels)
    val defaulted: ErrorOr[Seq[PipelinesApiAttachedDisk]] = addDefault(sequenced)
    defaulted
  }

  private def validateLocalDisk(disk: String): ErrorOr[PipelinesApiAttachedDisk] = {
    PipelinesApiAttachedDisk.parse(disk) match {
      case scala.util.Success(attachedDisk) => attachedDisk.validNel
      case scala.util.Failure(ex) => ex.getMessage.invalidNel
    }
  }

  private def sequenceNels(nels: Seq[ErrorOr[PipelinesApiAttachedDisk]]): ErrorOr[Seq[PipelinesApiAttachedDisk]] = {
    val emptyDiskNel: ErrorOr[Vector[PipelinesApiAttachedDisk]] = Vector.empty[PipelinesApiAttachedDisk].validNel
    val disksNel: ErrorOr[Vector[PipelinesApiAttachedDisk]] = nels.foldLeft(emptyDiskNel) {
      (acc, v) => (acc, v) mapN { (a, v) => a :+ v }
    }
    disksNel
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
