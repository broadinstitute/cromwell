package cromwell.backend.impl.jes

import cromwell.backend.impl.jes.io.{JesWorkingDisk, JesAttachedDisk}
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation.RuntimeAttributesValidation._
import cromwell.backend.validation._
import cromwell.core._
import lenthall.exception.MessageAggregation
import wdl4s.types.{WdlStringType, WdlIntegerType}
import wdl4s.values.{WdlInteger, WdlArray, WdlString, WdlValue}

import scalaz._
import Scalaz._

case class JesRuntimeAttributes(cpu: Int,
                                zones: Vector[String],
                                preemptible: Int,
                                bootDiskSize: Int,
                                memory: MemorySize,
                                disks: Seq[JesAttachedDisk],
                                dockerImage: Option[String],
                                failOnStderr: Boolean,
                                continueOnReturnCode: ContinueOnReturnCode)

object JesRuntimeAttributes {
  val CpuDefaultValue = 1
  val ZonesKey = "zones"
  val ZoneDefaultValue = Seq("us-central1-a")
  val PreemptibleKey = "preemptible"
  val PreemptibleDefaultValue = 0
  val BootDiskSizeKey = "bootDiskSizeGb"
  val BootDiskSizeDefaultValue = 10
  val MemoryDefaultValue = "2 GB"
  val DisksKey = "disks"
  val DisksDefaultValue = s"${JesWorkingDisk.Name} 10 SSD"
  val DefaultJesWorkingDisk = JesAttachedDisk.parse(DisksDefaultValue).get
  val FailOnStderrDefaultValue = false
  val ContinueOnRcDefaultValue = 0

  def apply(attrs: Map[String, WdlValue]): JesRuntimeAttributes = {
    val cpu = validateCpu(attrs.get(Cpu), CpuDefaultValue.successNel)
    val zones = validateZone(attrs.get(ZonesKey))
    val preemptible = validatePreemptible(attrs.get(PreemptibleKey))
    val bootDiskSize = validateBootDisk(attrs.get(BootDiskSizeKey))
    val memory = validateMemory(attrs.get(Memory), parseMemoryString(WdlString(MemoryDefaultValue)))
    val disks = validateLocalDisks(attrs.get(DisksKey))
    val docker = validateDocker(attrs.get(Docker), "Failed to get Docker mandatory key from runtime attributes".failureNel)
    val failOnStderr = validateFailOnStderr(attrs.get(FailOnStderr), FailOnStderrDefaultValue.successNel)
    val continueOnReturnCode = validateContinueOnReturnCode(attrs.get(ContinueOnReturnCode),
      ContinueOnReturnCodeSet(Set(ContinueOnRcDefaultValue)).successNel)
    (cpu |@| zones |@| preemptible |@| bootDiskSize |@| memory |@| disks |@| docker |@| failOnStderr |@| continueOnReturnCode) {
      new JesRuntimeAttributes(_, _, _, _, _, _, _, _, _)
    } match {
      case Success(x) => x
      case Failure(nel) => throw new RuntimeException with MessageAggregation {
        override def exceptionContext: String = "Runtime attribute validation failed"
        override def errorMessages: Traversable[String] = nel.list
      }
    }
  }

  private def validateZone(zoneValue: Option[WdlValue]): ErrorOr[Vector[String]] = {
    zoneValue match {
      case Some(WdlString(s)) => s.split("\\s+").toVector.successNel
      case Some(WdlArray(wdlType, value)) if wdlType.memberType == WdlStringType =>
        value.map(_.valueString).toVector.successNel
      case Some(_) => s"Expecting $ZonesKey runtime attribute to be either a whitespace separated String or an Array[String]".failureNel
      case None => ZoneDefaultValue.toVector.successNel
    }
  }

  private def validatePreemptible(preemptible: Option[WdlValue]): ErrorOr[Int] = {
    val preemptibleValidation = preemptible.map(validateInt).getOrElse(PreemptibleDefaultValue.successNel)
    if (preemptibleValidation.isFailure) {
      s"Expecting $PreemptibleKey runtime attribute to be an Integer".failureNel
    }
    else {
      preemptibleValidation
    }
  }

  private def validateBootDisk(diskSize: Option[WdlValue]): ErrorOr[Int] = diskSize match {
    case Some(x) if WdlIntegerType.isCoerceableFrom(x.wdlType) =>
      WdlIntegerType.coerceRawValue(x) match {
        case scala.util.Success(x: WdlInteger) => x.value.intValue.successNel
        case scala.util.Success(unhandled) => s"Coercion was expected to create an Integer but instead got $unhandled".failureNel
        case scala.util.Failure(t) => s"Expecting $BootDiskSizeKey runtime attribute to be an Integer".failureNel
      }
    case None => BootDiskSizeDefaultValue.successNel
  }

  private def validateLocalDisks(value: Option[WdlValue]): ErrorOr[Seq[JesAttachedDisk]] = {
    val nels = value match {
      case Some(WdlString(s)) => s.split(",\\s*").toSeq.map(validateLocalDisk)
      case Some(WdlArray(wdlType, seq)) if wdlType.memberType == WdlStringType =>
        seq.map(_.valueString).map(validateLocalDisk)
      case Some(_) =>
        Seq(s"Expecting $DisksKey runtime attribute to be a comma separated String or Array[String]".failureNel[JesAttachedDisk])
      case None => Seq(DefaultJesWorkingDisk).map(_.successNel)
    }

    val emptyDiskNel = Vector.empty[JesAttachedDisk].successNel[String]
    val disksNel = nels.foldLeft(emptyDiskNel)((acc, v) => (acc |@| v) { (a, v) => a :+ v })

    disksNel map {
      case disks if disks.exists(_.name == JesWorkingDisk.Name) => disks
      case disks => disks :+ DefaultJesWorkingDisk
    }
  }

  private def validateLocalDisk(disk: String): ErrorOr[JesAttachedDisk] = {
    JesAttachedDisk.parse(disk) match {
      case scala.util.Success(localDisk) => localDisk.successNel
      case scala.util.Failure(ex) => ex.getMessage.failureNel
    }
  }

}
