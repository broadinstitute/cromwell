package cromwell.engine.backend.runtimeattributes

import com.google.api.services.genomics.model.Disk
import cromwell.engine.backend.BackendType
import cromwell.engine.backend.runtimeattributes.RuntimeKey._
import cromwell.engine.workflow.WorkflowOptions
import org.slf4j.LoggerFactory
import wdl4s._
import wdl4s.parser.MemorySize

import scala.language.postfixOps
import scalaz.Scalaz._
import scalaz._

case class CromwellRuntimeAttributes(docker: Option[String],
                                     zones: Seq[String],
                                     failOnStderr: Boolean,
                                     continueOnReturnCode: ContinueOnReturnCode,
                                     cpu: Long,
                                     preemptible: Int,
                                     disks: Seq[Disk],
                                     memoryGB: Double)

object CromwellRuntimeAttributes {
  private val log = LoggerFactory.getLogger("RuntimeAttributes")

  def apply(wdlRuntimeAttributes: RuntimeAttributes, backendType: BackendType): CromwellRuntimeAttributes =
    fromAttributeMap(AttributeMap(wdlRuntimeAttributes.attrs), backendType)

  def apply(wdlRuntimeAttributes: RuntimeAttributes, workflowOptions: WorkflowOptions, backendType: BackendType) = {
    val attributesMap = AttributeMap(wdlRuntimeAttributes.attrs)

    val withDefaults = new AttributeMapTrait {
      override def keys = attributesMap.keys ++ workflowOptions.getDefaultRuntimeOptionKeys
      override def unsupportedKeys(bt: BackendType): Seq[String] = attributesMap.unsupportedKeys(bt)
      override def getSeq(key: RuntimeKey): Option[Seq[String]] = attributesMap.getSeq(key)
      override def get(key: RuntimeKey): Option[String] = attributesMap.get(key) match {
        case Some(x) => Some(x)
        case None => workflowOptions.getDefaultRuntimeOption(key.key).toOption
      }
    }

    fromAttributeMap(withDefaults, backendType)
  }

  private def fromAttributeMap(attributeMap: AttributeMapTrait, backendType: BackendType): CromwellRuntimeAttributes = {
    /**
      *  Warn if keys are found that are unsupported on this backend.  This is not necessarily an error, at this point
      *  the keys are known to be supported on at least one backend other than this.
      */
    attributeMap.unsupportedKeys(backendType).foreach(unsupportedKeyMessage => log.info(unsupportedKeyMessage))

    val attributeMapNel = validateRequiredKeys(attributeMap, backendType)
    val runtimeAttributeNel = validateRuntimeAttributes(attributeMap)
    val validatedRuntimeAttributes = (attributeMapNel |@| runtimeAttributeNel) { (_, r) => r }

    validatedRuntimeAttributes match {
      case Success(r) => r
      case Failure(f) =>
        throw new IllegalArgumentException() with ThrowableWithErrors {
          val message = "RuntimeAttribute is not valid."
          val errors = f
        }
    }
  }

  val LocalDiskName = "local-disk"
  val LocalizationDisk = LocalDisk(LocalDiskName, DiskType.SSD, Option(10)).toDisk

  /** Fallback values if values for these keys are not specified in a task "runtime" stanza. */
  object Defaults {
    val Cpu = 1L
    val Disk = Vector(LocalizationDisk)
    val FailOnStderr = false
    val ContinueOnReturnCode = ContinueOnReturnCodeSet(Set(0))
    val Memory = 2.0
    val Preemptible = 0
    val Zones = Vector("us-central1-a")
  }

  private def validateRuntimeAttributes(attributeMap: AttributeMapTrait): ValidationNel[String, CromwellRuntimeAttributes] = {
    val docker = attributeMap.get(DOCKER)
    val zones = attributeMap.get(ZONES) map { _.split("\\s+").toVector } getOrElse Defaults.Zones
    val failOnStderr = validateFailOnStderr(attributeMap.get(FAIL_ON_STDERR))
    val continueOnReturnCode = validateContinueOnReturnCode(attributeMap.getSeq(CONTINUE_ON_RETURN_CODE))
    val cpu = validateCpu(attributeMap.get(CPU))
    val disks = validateLocalDisks(attributeMap.get(DISKS))
    val preemptible = validatePreemptible(attributeMap.get(PREEMPTIBLE))
    val memory = validateMemory(attributeMap.get(MEMORY))

    (failOnStderr |@| continueOnReturnCode |@| cpu |@| preemptible |@| disks |@| memory) {
      CromwellRuntimeAttributes(docker, zones, _, _, _, _, _, _)
    }
  }

  private def validateCpu(value: Option[String]): ValidationNel[String, Long] = {
    value map validateLong getOrElse Defaults.Cpu.successNel
  }

  private def validatePreemptible(value: Option[String]): ValidationNel[String, Int] = {
    value map validateInt getOrElse Defaults.Preemptible.successNel
  }

  private def validateFailOnStderr(value: Option[String]): ValidationNel[String, Boolean] = {
    value map validateBoolean getOrElse Defaults.FailOnStderr.successNel
  }

  private def validateContinueOnReturnCode(value: Option[Seq[String]]): ValidationNel[String, ContinueOnReturnCode] = {
    value.map(_.map(_.toLowerCase)) match {
      case Some(Seq("true")) => ContinueOnReturnCodeFlag(true).successNel
      case Some(Seq("false")) => ContinueOnReturnCodeFlag(false).successNel
      case Some(seq) =>
        val nels = seq map validateContinueOnReturnCode
        val defaultReturnCodeNel = Set.empty[Int].successNel[String]
        nels.foldLeft(defaultReturnCodeNel)((acc, v) => (acc |@| v) { (a, v) => a + v }) map ContinueOnReturnCodeSet
      case None => Defaults.ContinueOnReturnCode.successNel
    }
  }

  private def validateContinueOnReturnCode(returnCode: String): ValidationNel[String, Int] = validateInt(returnCode)

  private def validateLocalDisks(value: Option[String]): ValidationNel[String, Seq[Disk]] = {
    value match {
      case Some(v) =>
        val nels = v.split(",\\s*") map validateLocalDisk
        val emptyDiskNel = Vector.empty[Disk].successNel[String]
        val disksNel = nels.foldLeft(emptyDiskNel)((acc, v) => (acc |@| v) { (a, v) => a :+ v })
        disksNel map {
          case disks if disks.exists(_.getName == LocalDiskName) => disks
          case disks => disks :+ LocalizationDisk
        }
      case None => Defaults.Disk.successNel
    }
  }

  private def validateLocalDisk(disk: String): ValidationNel[String, Disk] = {
    disk.split("\\s+") match {
      case Array(name, DiskType.LOCAL.diskTypeName) =>
        LocalDisk(name, DiskType.LOCAL).toDisk.successNel[String]
      case Array(name, sizeGb, diskType) if diskType != DiskType.LOCAL.diskTypeName =>
        (validateLong(sizeGb) |@| validateDiskType(diskType)) { (s, dt) => LocalDisk(name, dt, Option(s)).toDisk }
      case _ => s"'$disk' should be in form 'NAME SIZE TYPE', with SIZE blank for LOCAL, otherwise SIZE in GB".failureNel
    }
  }

  private def validateDiskType(diskTypeName: String): ValidationNel[String, DiskType] = {
    DiskType.values().find(_.diskTypeName == diskTypeName) match {
      case Some(diskType) => diskType.successNel[String]
      case None =>
        val diskTypeNames = DiskType.values.map(_.diskTypeName).mkString(", ")
        s"Disk TYPE $diskTypeName should be one of $diskTypeNames".failureNel
    }
  }

  private def validateMemory(value: Option[String]): ValidationNel[String, Double] = {
    value map validateMemStringInGb getOrElse Defaults.Memory.successNel
  }

  private def validateMemStringInGb(mem: String): ValidationNel[String, Double] = {
    mem match {
      case MemoryAndUnitPattern(amountString, unitString) =>
        val amount = validateLong(amountString)
        val unit = validateMemoryUnit(unitString)
        (amount |@| unit) { (a, u) => MemorySize.GB.fromBytes(u.toBytes(a)) }
      case _ => s"$mem should be of the form X Unit where X is a number, e.g. 8 GB".failureNel
    }
  }

  private def validateMemoryUnit(unit: String): ValidationNel[String, MemorySize] = {
    MemorySize.values find { _.suffixes.contains(unit) } match {
      case Some(s) => s.successNel
      case None => s"$unit is an invalid memory unit".failureNel
    }
  }

  private def validateRequiredKeys(attributeMap: AttributeMapTrait,
                                   backendType: BackendType): ValidationNel[String, Unit] = {
    val requiredKeys = for {
      key <- RuntimeKey.values().toSet
      if key.isMandatory(backendType)
    } yield key.key

    val missingRequiredKeys = requiredKeys -- attributeMap.keys
    if (missingRequiredKeys.isEmpty) ().successNel
    else {
      val missingKeyString = missingRequiredKeys.toSeq.sorted.mkString(", ")
      s"Missing required keys in runtime configuration for backend '$backendType': $missingKeyString".failureNel
    }
  }

  private def unknownKeys(attributeMap: AttributeMapTrait,
                                  backendType: BackendType): Set[String] = {
    val knownKeys = RuntimeKey.values map { _.key }
    // Finding keys unknown to any backend is an error.
    attributeMap.keys -- knownKeys
  }

  private def validateBoolean(value: String): ValidationNel[String, Boolean] = {
    try {
      value.toBoolean.successNel
    } catch {
      case _: IllegalArgumentException => s"$value not convertible to a Boolean".failureNel[Boolean]
    }
  }

  private def validateInt(value: String): ValidationNel[String, Int] = {
    try {
      value.toInt.successNel
    } catch {
      case _: IllegalArgumentException => s"$value not convertible to an Int".failureNel[Int]
    }
  }

  private def validateLong(value: String): ValidationNel[String, Long] = {
    try {
      value.toLong.successNel
    } catch {
      case _: IllegalArgumentException => s"$value not convertible to a Long".failureNel[Long]
    }
  }

  private case class LocalDisk(name: String, diskType: DiskType, sizeGbOption: Option[Long] = None) {
    def toDisk: Disk = {
      val disk = new Disk().setName(name).setType(diskType.googleTypeName).setAutoDelete(true)
      // Even though GCE ignores the value for local disks, JES requires we set the disk size anyway.
      disk.setSizeGb(long2Long(sizeGbOption getOrElse 10L))
      disk
    }
  }

  implicit class EnhancedBackendType(val backendType: BackendType) extends AnyVal {
    def supportedKeys: Set[RuntimeKey] = for {
      key <- RuntimeKey.values().toSet
      if key.supports(backendType)
    } yield key
  }

  implicit class EnhancedAttributeMap(val attrs: Map[String, String]) extends AnyVal {
    def unsupportedKeys(backendType: BackendType): Seq[String] = {
      val supportedKeys = backendType.supportedKeys map { _.key }
      val unsupportedKeys = attrs.keySet -- supportedKeys

      if (unsupportedKeys.isEmpty) Vector.empty
      else Vector(s"Found unsupported keys for backend '$backendType': " + unsupportedKeys.toSeq.sorted.mkString(", "))
    }
  }

  val DefaultReturnCodeNel = Set.empty[Int].successNel[String]
  val MemoryAndUnitPattern = """(\d+)\s*(\w+)""".r
}
