package cromwell.binding

import com.google.api.services.genomics.model.Disk
import cromwell.binding.AstTools.{AstNodeName, EnhancedAstNode}
import cromwell.binding.expression.NoFunctions
import cromwell.binding.values._
import cromwell.parser.RuntimeKey._
import cromwell.parser.WdlParser.{Ast, AstList}
import cromwell.parser.{BackendType, MemorySize, RuntimeKey}
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.language.postfixOps
import scalaz.Scalaz._
import scalaz._

case class RuntimeAttributes(docker: Option[String],
              defaultZones: Seq[String],
              failOnStderr: Boolean,
              continueOnReturnCode: ContinueOnReturnCode,
              cpu: Long,
              preemptible: Boolean,
              defaultDisks: Seq[Disk],
              memoryGB: Double)

object RuntimeAttributes {
  private val log = LoggerFactory.getLogger("RuntimeAttributes")

  def apply(ast: Ast, backendType: BackendType): RuntimeAttributes = {
    val attributeMap = ast.toAttributes

    /**
     *  Warn if keys are found that are unsupported on this backend.  This is not necessarily an error, at this point
     *  the keys are known to be supported on at least one backend other than this.
     */
    attributeMap.unsupportedKeys(backendType) foreach log.warn

    val attributeMapNel = validateAttributeMap(attributeMap, backendType)
    val runtimeAttributeNel = validateRuntimeAttributes(attributeMap)
    val validatedRuntimeAttributes = (attributeMapNel |@| runtimeAttributeNel) { (_, r) => r }

    validatedRuntimeAttributes match {
      case Success(r) => r
      case Failure(f) =>
        val errorMessages = f.list.mkString(", ")
        throw new IllegalArgumentException(s"RuntimeAttribute is not valid: Errors: $errorMessages")
    }
  }

  val LocalDiskName = "local-disk"

  val LocalizationDisk = LocalDisk(LocalDiskName, DiskType.LOCAL).toDisk
  
  /** Fallback values if values for these keys are not specified in a task "runtime" stanza. */
  object Defaults {
    val Cpu = 1L
    val Disk = Vector(LocalizationDisk)
    val FailOnStderr = false
    val ContinueOnReturnCode = ContinueOnReturnCodeSet(Set(0))
    val Memory = 2.0
    val Preemptible = false
    val Zones = Vector("us-central1-a")
  }

  private def validateRuntimeAttributes(attributeMap: AttributeMap): ValidationNel[String, RuntimeAttributes] = {
    val docker = attributeMap.get(DOCKER)
    val zones = attributeMap.get(DEFAULT_ZONES) map { _.split("\\s+").toVector } getOrElse Defaults.Zones
    val failOnStderr = validateFailOnStderr(attributeMap.get(FAIL_ON_STDERR))
    val continueOnReturnCode = validateContinueOnReturnCode(attributeMap.getSeq(CONTINUE_ON_RETURN_CODE))
    val cpu = validateCpu(attributeMap.get(CPU))
    val preemptible = validatePreemptible(attributeMap.get(PREEMPTIBLE))
    val disks = validateLocalDisks(attributeMap.get(DEFAULT_DISKS))
    val memory = validateMemory(attributeMap.get(MEMORY))

    (failOnStderr |@| continueOnReturnCode |@| cpu |@| preemptible |@| disks |@| memory) {
      RuntimeAttributes(docker, zones, _, _, _, _, _, _)
    }
  }

  private def validateCpu(value: Option[String]): ValidationNel[String, Long] = {
    value map validateLong getOrElse Defaults.Cpu.successNel
  }

  private def validatePreemptible(value: Option[String]): ValidationNel[String, Boolean] = {
    value map validateBoolean getOrElse Defaults.Preemptible.successNel
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

  private def validateContinueOnReturnCode(returnCode: String): ValidationNel[String, Int] = {
    validateInt(returnCode)
  }

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
    val memoryPattern = """(\d+)\s*(\w+)""".r
    mem match {
      case memoryPattern(amountString, unitString) =>
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

  private def validateAttributeMap(attributeMap: AttributeMap,
                                   backendType: BackendType): ValidationNel[String, Unit] = {
    val requiredKeysNel = validateRequiredKeys(attributeMap, backendType)
    val unknownKeysNel = validateUnknownKeys(attributeMap, backendType)
    /*
      Combining either the two FailureNELs and their respective error message or "combining" the two
      SuccessNELs and their units. The latter is really a side effect to carry the success forward
     */
    (requiredKeysNel |@| unknownKeysNel) { (_, _) }
  }

  private def validateRequiredKeys(attributeMap: AttributeMap,
                                   backendType: BackendType): ValidationNel[String, Unit] = {
    val requiredKeys = for {
      key <- RuntimeKey.values().toSet
      if key.isMandatory(backendType)
    } yield key.key

    val missingRequiredKeys = requiredKeys -- attributeMap.attrs.keySet
    if (missingRequiredKeys.isEmpty) ().successNel
    else {
      val missingKeyString = missingRequiredKeys.toSeq.sorted.mkString(", ")
      s"Missing required keys in runtime configuration for backend '$backendType': $missingKeyString".failureNel
    }
  }

  private def validateUnknownKeys(attributeMap: AttributeMap,
                                  backendType: BackendType): ValidationNel[String, Unit] = {
    val knownKeys = RuntimeKey.values map { _.key }
    // Finding keys unknown to any backend is an error.
    val unknownKeys = attributeMap.attrs.keySet -- knownKeys

    if (unknownKeys.isEmpty) ().successNel
    else {
      val unknownKeyString = unknownKeys.toSeq.sorted.mkString(", ")
      s"Unknown keys found in runtime configuration: $unknownKeyString".failureNel
    }
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

  private def processRuntimeAttribute(ast: Ast): (String, Seq[String]) = {
    val key = ast.getAttribute("key").sourceString
    val seq = Option(ast.getAttribute("value")) map { valAttr =>
      // TODO: NoFunctions should be converted to a case object or a constant
      WdlExpression.evaluate(valAttr, NoLookup, new NoFunctions) match {
        case scala.util.Success(wdlPrimitive: WdlPrimitive) => Seq(wdlPrimitive.valueString)
        case scala.util.Success(wdlArray: WdlArray) => wdlArray.value.map(_.valueString)
        case scala.util.Success(wdlValue: WdlValue) =>
          throw new IllegalArgumentException(
            s"WdlType not supported for $key, ${wdlValue.wdlType}: ${wdlValue.valueString}")
        case null =>
          throw new IllegalArgumentException(s"value was null: $key}")
        case scala.util.Failure(f) => throw f
      }
    } getOrElse Seq.empty
    (key, seq)
  }

  private def processRuntimeAttributes(astList: AstList): Map[String, Seq[String]] = {
    astList.asScala.toVector map { a => processRuntimeAttribute(a.asInstanceOf[Ast]) } toMap
  }

  private case class LocalDisk(name: String, diskType: DiskType, sizeGbOption: Option[Long] = None) {
    def toDisk: Disk = {
      val disk = new Disk().setName(name).setType(diskType.googleTypeName).setAutoDelete(true)
      // Even though GCE ignores the value for local disks, JES requires we set the disk size anyway.
      disk.setSizeGb(long2Long(sizeGbOption getOrElse 10L))
      disk
    }
  }

  implicit class EnhancedAst(val ast: Ast) extends AnyVal {
    def toAttributes: AttributeMap = {
      val asts = ast.findAsts(AstNodeName.Runtime)
      if (asts.size > 1) throw new UnsupportedOperationException("Only one runtime block may be defined per task")
      val astList = asts.headOption map { _.getAttribute("map").asInstanceOf[AstList] }
      AttributeMap(astList map processRuntimeAttributes getOrElse Map.empty[String, Seq[String]])
    }
  }

  val MemoryAndUnitPattern = "(\\d+)\\s*(\\w+)".r
}
