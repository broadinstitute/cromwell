package cromwell.backend.impl.jes

import cats.data.NonEmptyList
import cats.data.Validated._
import cats.syntax.cartesian._
import cats.syntax.validated._
import cromwell.backend.MemorySize
import cromwell.backend.impl.jes.io.{JesAttachedDisk, JesWorkingDisk}
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation.{BooleanRuntimeAttributesValidation, _}
import lenthall.validation.ErrorOr._
import wdl4s.types._
import wdl4s.values._

case class JesRuntimeAttributes(cpu: Int,
                                zones: Vector[String],
                                preemptible: Int,
                                bootDiskSize: Int,
                                memory: MemorySize,
                                disks: Seq[JesAttachedDisk],
                                dockerImage: String,
                                failOnStderr: Boolean,
                                continueOnReturnCode: ContinueOnReturnCode,
                                noAddress: Boolean)

object JesRuntimeAttributes {
  private val MemoryDefaultValue = "2 GB"

  val ZonesKey = "zones"

  val PreemptibleKey = "preemptible"
  private val PreemptibleDefaultValue = 0

  val BootDiskSizeKey = "bootDiskSizeGb"
  private val BootDiskSizeDefaultValue = 10

  val NoAddressKey = "noAddress"
  private val NoAddressDefaultValue = false

  val DisksKey = "disks"
  private val DisksDefaultValue = s"${JesWorkingDisk.Name} 10 SSD"

  private val cpuValidation: RuntimeAttributesValidation[Int] = CpuValidation.default

  private val disksValidation: RuntimeAttributesValidation[Seq[JesAttachedDisk]] =
    DisksValidation.withDefault(WdlString(DisksDefaultValue))

  private def zonesValidation(defaultZones: NonEmptyList[String]): RuntimeAttributesValidation[Vector[String]] =
    ZonesValidation.withDefault(WdlString(defaultZones.toList.mkString(" ")))

  private val preemptibleValidation: RuntimeAttributesValidation[Int] =
    new IntRuntimeAttributesValidation(JesRuntimeAttributes.PreemptibleKey)
      .withDefault(WdlInteger(PreemptibleDefaultValue))

  private val memoryValidation: RuntimeAttributesValidation[MemorySize] =
    MemoryValidation.withDefaultMemory(MemorySize.parse(MemoryDefaultValue).get)

  private val bootDiskSizeValidation: RuntimeAttributesValidation[Int] =
    new IntRuntimeAttributesValidation(JesRuntimeAttributes.BootDiskSizeKey)
      .withDefault(WdlInteger(BootDiskSizeDefaultValue))

  private val noAddressValidation: RuntimeAttributesValidation[Boolean] =
    new BooleanRuntimeAttributesValidation(JesRuntimeAttributes.NoAddressKey)
      .withDefault(WdlBoolean(NoAddressDefaultValue))

  private val dockerValidation: RuntimeAttributesValidation[String] = DockerValidation.instance

  def runtimeAttributesBuilder(jesConfiguration: JesConfiguration): StandardValidatedRuntimeAttributesBuilder =
    StandardValidatedRuntimeAttributesBuilder.default.withValidation(
      cpuValidation,
      disksValidation,
      zonesValidation(jesConfiguration.defaultZones),
      preemptibleValidation,
      memoryValidation,
      bootDiskSizeValidation,
      noAddressValidation,
      dockerValidation
    )

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes): JesRuntimeAttributes = {
    val cpu: Int = RuntimeAttributesValidation.extract(cpuValidation, validatedRuntimeAttributes)
    val zones: Vector[String] = RuntimeAttributesValidation.extract(ZonesValidation, validatedRuntimeAttributes)
    val preemptible: Int = RuntimeAttributesValidation.extract(preemptibleValidation, validatedRuntimeAttributes)
    val bootDiskSize: Int = RuntimeAttributesValidation.extract(bootDiskSizeValidation, validatedRuntimeAttributes)
    val memory: MemorySize = RuntimeAttributesValidation.extract(memoryValidation, validatedRuntimeAttributes)
    val disks: Seq[JesAttachedDisk] = RuntimeAttributesValidation.extract(disksValidation, validatedRuntimeAttributes)
    val docker: String = RuntimeAttributesValidation.extract(dockerValidation, validatedRuntimeAttributes)
    val failOnStderr: Boolean =
      RuntimeAttributesValidation.extract(FailOnStderrValidation.default, validatedRuntimeAttributes)
    val continueOnReturnCode: ContinueOnReturnCode =
      RuntimeAttributesValidation.extract(ContinueOnReturnCodeValidation.default, validatedRuntimeAttributes)
    val noAddress: Boolean = RuntimeAttributesValidation.extract(noAddressValidation, validatedRuntimeAttributes)

    new JesRuntimeAttributes(
      cpu,
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
  override def key: String = JesRuntimeAttributes.ZonesKey

  override def coercion: Traversable[WdlType] = Set(WdlStringType, WdlArrayType(WdlStringType))

  override protected def validateValue: PartialFunction[WdlValue, ErrorOr[Vector[String]]] = {
    case WdlString(s) => s.split("\\s+").toVector.validNel
    case WdlArray(wdlType, value) if wdlType.memberType == WdlStringType =>
      value.map(_.valueString).toVector.validNel
  }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be either a whitespace separated String or an Array[String]"
}

object DisksValidation extends RuntimeAttributesValidation[Seq[JesAttachedDisk]] {
  override def key: String = JesRuntimeAttributes.DisksKey

  override def coercion: Traversable[WdlType] = Set(WdlStringType, WdlArrayType(WdlStringType))

  override protected def validateValue: PartialFunction[WdlValue, ErrorOr[Seq[JesAttachedDisk]]] = {
    case WdlString(value) => validateLocalDisks(value.split(",\\s*").toSeq)
    case WdlArray(wdlType, values) if wdlType.memberType == WdlStringType =>
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
    val emptyDiskNel = Vector.empty[JesAttachedDisk].validNel[String]
    val disksNel: ErrorOr[Vector[JesAttachedDisk]] = nels.foldLeft(emptyDiskNel) {
      (acc, v) => (acc |@| v) map { (a, v) => a :+ v }
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
