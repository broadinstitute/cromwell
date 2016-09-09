package cromwell.backend.impl.jes

import cromwell.backend.MemorySize
import cromwell.backend.impl.jes.io.{JesWorkingDisk, JesAttachedDisk}
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation.RuntimeAttributesValidation._
import cromwell.backend.validation._
import cromwell.core._
import lenthall.exception.MessageAggregation
import org.slf4j.Logger
import wdl4s.types._
import wdl4s.values._
import cromwell.backend.validation.RuntimeAttributesDefault._
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
                                continueOnReturnCode: ContinueOnReturnCode,
                                noAddress: Boolean) {
  import JesRuntimeAttributes._

  lazy val asMap = Map[String, Any](
    CpuKey -> cpu.toString,
    ZonesKey -> zones.mkString(","),
    PreemptibleKey -> preemptible.toString,
    BootDiskSizeKey -> bootDiskSize.toString,
    MemoryKey -> memory.toString,
    DisksKey -> disks.mkString(","),
    DockerKey -> dockerImage.get,
    FailOnStderrKey -> failOnStderr.toString,
    ContinueOnReturnCodeKey -> continueOnReturnCode
  )
}

object JesRuntimeAttributes {
  private val CpuDefaultValue = 1
  private val ContinueOnReturnCodeDefaultValue = 0
  private val MemoryDefaultValue = "2 GB"

  val ZonesKey = "zones"
  private val ZoneDefaultValue = "us-central1-a"

  val PreemptibleKey = "preemptible"
  private val PreemptibleDefaultValue = 0

  val BootDiskSizeKey = "bootDiskSizeGb"
  private val BootDiskSizeDefaultValue = 10

  val NoAddressKey = "noAddress"
  private val NoAddressDefaultValue = false

  val DisksKey = "disks"
  private val DisksDefaultValue = s"${JesWorkingDisk.Name} 10 SSD"

  val staticDefaults = Map(
    CpuKey -> WdlInteger(CpuDefaultValue),
    DisksKey -> WdlString(DisksDefaultValue),
    ZonesKey -> WdlString(ZoneDefaultValue),
    ContinueOnReturnCodeKey -> WdlInteger(ContinueOnReturnCodeDefaultValue),
    FailOnStderrKey -> WdlBoolean.False,
    PreemptibleKey -> WdlInteger(PreemptibleDefaultValue),
    MemoryKey -> WdlString(MemoryDefaultValue),
    BootDiskSizeKey -> WdlInteger(BootDiskSizeDefaultValue),
    NoAddressKey -> WdlBoolean(NoAddressDefaultValue)
  )

  private[jes] val coercionMap: Map[String, Set[WdlType]] = Map(
    CpuKey -> Set(WdlIntegerType),
    DisksKey -> Set(WdlStringType, WdlArrayType(WdlStringType)),
    ZonesKey -> Set(WdlStringType, WdlArrayType(WdlStringType)),
    ContinueOnReturnCodeKey -> ContinueOnReturnCode.validWdlTypes,
    FailOnStderrKey -> Set(WdlBooleanType),
    PreemptibleKey -> Set(WdlIntegerType),
    MemoryKey -> Set(WdlStringType),
    BootDiskSizeKey -> Set(WdlIntegerType),
    NoAddressKey -> Set(WdlBooleanType),
    DockerKey -> Set(WdlStringType)
  )

  def apply(attrs: Map[String, WdlValue], logger: Logger): JesRuntimeAttributes = {
    warnUnrecognized(attrs.keySet, coercionMap.keySet, logger)

    val cpu = validateCpu(attrs.get(CpuKey), noValueFoundFor(CpuKey))
    val memory = validateMemory(attrs.get(MemoryKey), noValueFoundFor(MemoryKey))
    val docker = validateDocker(attrs.get(DockerKey), noValueFoundFor(DockerKey))
    val failOnStderr = validateFailOnStderr(attrs.get(FailOnStderrKey), noValueFoundFor(FailOnStderrKey))
    val continueOnReturnCode = validateContinueOnReturnCode(attrs.get(ContinueOnReturnCodeKey), noValueFoundFor(ContinueOnReturnCodeKey))

    val zones = validateZone(attrs(ZonesKey))
    val preemptible = validatePreemptible(attrs(PreemptibleKey))
    val noAddress = validateNoAddress(attrs(NoAddressKey))
    val bootDiskSize = validateBootDisk(attrs(BootDiskSizeKey))
    val disks = validateLocalDisks(attrs(DisksKey))
    (cpu |@| zones |@| preemptible |@| bootDiskSize |@| memory |@| disks |@| docker |@| failOnStderr |@| continueOnReturnCode |@| noAddress) {
      new JesRuntimeAttributes(_, _, _, _, _, _, _, _, _, _)
    } match {
      case Success(x) => x
      case Failure(nel) => throw new RuntimeException with MessageAggregation {
        override def exceptionContext: String = "Runtime attribute validation failed"
        override def errorMessages: Traversable[String] = nel.list.toList
      }
    }
  }

  private def validateZone(zoneValue: WdlValue): ErrorOr[Vector[String]] = {
    zoneValue match {
      case WdlString(s) => s.split("\\s+").toVector.successNel
      case WdlArray(wdlType, value) if wdlType.memberType == WdlStringType =>
        value.map(_.valueString).toVector.successNel
      case _ => s"Expecting $ZonesKey runtime attribute to be either a whitespace separated String or an Array[String]".failureNel
    }
  }

  private def contextualizeFailure[T](validation: ErrorOr[T], key: String): ErrorOr[T] = {
    validation.leftMap[String](errors => s"Failed to validate $key runtime attribute: " + errors.toList.mkString(",")).toValidationNel
  }

  private def validatePreemptible(preemptible: WdlValue): ErrorOr[Int] = {
    contextualizeFailure(validateInt(preemptible), PreemptibleKey)
  }

  private def validateNoAddress(noAddress: WdlValue): ErrorOr[Boolean] = {
    contextualizeFailure(validateBoolean(noAddress), NoAddressKey)
  }

  private def validateBootDisk(diskSize: WdlValue): ErrorOr[Int] = diskSize match {
    case x if WdlIntegerType.isCoerceableFrom(x.wdlType) =>
      WdlIntegerType.coerceRawValue(x) match {
        case scala.util.Success(x: WdlInteger) => x.value.intValue.successNel
        case scala.util.Success(unhandled) => s"Coercion was expected to create an Integer but instead got $unhandled".failureNel
        case scala.util.Failure(t) => s"Expecting $BootDiskSizeKey runtime attribute to be an Integer".failureNel
      }
  }

  private def validateLocalDisks(value: WdlValue): ErrorOr[Seq[JesAttachedDisk]] = {
    val nels = value match {
      case WdlString(s) => s.split(",\\s*").toSeq.map(validateLocalDisk)
      case WdlArray(wdlType, seq) if wdlType.memberType == WdlStringType =>
        seq.map(_.valueString).map(validateLocalDisk)
      case _ =>
        Seq(s"Expecting $DisksKey runtime attribute to be a comma separated String or Array[String]".failureNel[JesAttachedDisk])
    }

    val emptyDiskNel = Vector.empty[JesAttachedDisk].successNel[String]
    val disksNel = nels.foldLeft(emptyDiskNel)((acc, v) => (acc |@| v) { (a, v) => a :+ v })

    disksNel map {
      case disks if disks.exists(_.name == JesWorkingDisk.Name) => disks
      case disks => disks :+ JesAttachedDisk.parse(DisksDefaultValue).get
    }
  }

  private def validateLocalDisk(disk: String): ErrorOr[JesAttachedDisk] = {
    JesAttachedDisk.parse(disk) match {
      case scala.util.Success(localDisk) => localDisk.successNel
      case scala.util.Failure(ex) => ex.getMessage.failureNel
    }
  }

}
