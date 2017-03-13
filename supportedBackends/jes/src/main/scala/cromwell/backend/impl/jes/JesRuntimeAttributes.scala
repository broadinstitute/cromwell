package cromwell.backend.impl.jes

import cats.data.NonEmptyList
import cats.syntax.validated._
import cats.data.Validated._
import cats.syntax.cartesian._
import com.typesafe.config.Config
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

  val ZonesKey = "zones"

  val PreemptibleKey = "preemptible"
  private val preemptibleValidationInstance = new IntRuntimeAttributesValidation(JesRuntimeAttributes.PreemptibleKey)

  val BootDiskSizeKey = "bootDiskSizeGb"
  private val bootDiskValidationInstance = new IntRuntimeAttributesValidation(JesRuntimeAttributes.BootDiskSizeKey)

  val NoAddressKey = "noAddress"
  private val noAddressValidationInstance = new BooleanRuntimeAttributesValidation(JesRuntimeAttributes.NoAddressKey)

  val DisksKey = "disks"

  private def cpuValidation(runtimeConfig: Config): RuntimeAttributesValidation[Int] =
    CpuValidation.default(runtimeConfig)

  private def disksValidation(runtimeConfig: Config): RuntimeAttributesValidation[Seq[JesAttachedDisk]] = DisksValidation
    .withDefault(DisksValidation.configDefaultWdlValue(runtimeConfig))

  private def zonesValidation(defaultZones: NonEmptyList[String]): RuntimeAttributesValidation[Vector[String]] = ZonesValidation
    .withDefault(Some(WdlString(defaultZones.toList.mkString(" "))))

  private def preemptibleValidation(runtimeConfig: Config): RuntimeAttributesValidation[Int] = preemptibleValidationInstance
    .withDefault(preemptibleValidationInstance.configDefaultWdlValue(runtimeConfig))

  private def memoryValidation(runtimeConfig: Config): RuntimeAttributesValidation[MemorySize] = {
    val memoryDefault = MemoryValidation.configDefaultString(runtimeConfig)
    MemoryValidation.withDefaultMemory(MemorySize.parse(memoryDefault), memoryDefault)
  }

  private def bootDiskSizeValidation(runtimeConfig: Config): RuntimeAttributesValidation[Int] = bootDiskValidationInstance
    .withDefault(bootDiskValidationInstance.configDefaultWdlValue(runtimeConfig))

  private def noAddressValidation(runtimeConfig: Config): RuntimeAttributesValidation[Boolean] = noAddressValidationInstance
    .withDefault(noAddressValidationInstance.configDefaultWdlValue(runtimeConfig))

  private val dockerValidation: RuntimeAttributesValidation[String] = DockerValidation.instance

  def runtimeAttributesBuilder(jesConfiguration: JesConfiguration): StandardValidatedRuntimeAttributesBuilder = {
    val runtimeConfig = jesConfiguration.runtimeAttributesConfig
    StandardValidatedRuntimeAttributesBuilder.default(runtimeConfig).withValidation(
      cpuValidation(runtimeConfig),
      disksValidation(runtimeConfig),
      zonesValidation(jesConfiguration.defaultZones),
      preemptibleValidation(runtimeConfig),
      memoryValidation(runtimeConfig),
      bootDiskSizeValidation(runtimeConfig),
      noAddressValidation(runtimeConfig),
      dockerValidation
    )
  }


  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes, runtimeAttrsConfig: Config): JesRuntimeAttributes = {
    val cpu: Int = RuntimeAttributesValidation.extract(cpuValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val zones: Vector[String] = RuntimeAttributesValidation.extract(ZonesValidation, validatedRuntimeAttributes)
    val preemptible: Int = RuntimeAttributesValidation.extract(preemptibleValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val bootDiskSize: Int = RuntimeAttributesValidation.extract(bootDiskSizeValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val memory: MemorySize = RuntimeAttributesValidation.extract(memoryValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val disks: Seq[JesAttachedDisk] = RuntimeAttributesValidation.extract(disksValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val docker: String = RuntimeAttributesValidation.extract(dockerValidation, validatedRuntimeAttributes)
    val failOnStderr: Boolean = RuntimeAttributesValidation.extract(FailOnStderrValidation.default(runtimeAttrsConfig), validatedRuntimeAttributes)
    val continueOnReturnCode: ContinueOnReturnCode = RuntimeAttributesValidation.extract(ContinueOnReturnCodeValidation.default(runtimeAttrsConfig), validatedRuntimeAttributes)
    val noAddress: Boolean = RuntimeAttributesValidation.extract(noAddressValidation(runtimeAttrsConfig), validatedRuntimeAttributes)

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
