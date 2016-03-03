package cromwell.engine.backend.runtimeattributes

import com.google.api.services.genomics.model.Disk
import cromwell.engine.ErrorOr
import cromwell.engine.backend.runtimeattributes.RuntimeKey._
import cromwell.engine.backend.{BackendCall, BackendType}
import cromwell.engine.workflow.WorkflowOptions
import cromwell.util.TryUtil
import org.slf4j.LoggerFactory
import wdl4s._
import wdl4s.parser.MemoryUnit
import wdl4s.types._
import wdl4s.values._

import scala.language.postfixOps
import scala.util.Try
import scalaz.Scalaz._
import scalaz._

case class CromwellRuntimeAttributes(attributes: Map[String, WdlValue],
                                     docker: Option[String],
                                     zones: Seq[String],
                                     failOnStderr: Boolean,
                                     continueOnReturnCode: ContinueOnReturnCode,
                                     cpu: Long,
                                     preemptible: Int,
                                     disks: Seq[Disk],
                                     memoryGB: Double)

object CromwellRuntimeAttributes {
  private val log = LoggerFactory.getLogger("RuntimeAttributes")

  def apply(wdlRuntimeAttributes: RuntimeAttributes, backendCall: BackendCall, workflowOptions: Option[WorkflowOptions]): CromwellRuntimeAttributes = {
    val supportedKeys = backendCall.backend.backendType.supportedKeys map { _.key }

    val attributes = for {
      attributesFromTask <- TryUtil.sequenceMap(wdlRuntimeAttributes.evaluate(backendCall.lookupFunction(backendCall.locallyQualifiedInputs), backendCall.engineFunctions))
      attributesWithDefaults <- Try(getAttributesWithDefaults(attributesFromTask, workflowOptions))
      _ <- validateKeys(attributesWithDefaults.keySet, backendCall.backend.backendType)
      supportedAttributes = attributesWithDefaults.filterKeys(k => supportedKeys.contains(k))
      validatedAttributes <- validateRuntimeAttributes(supportedAttributes)
    } yield validatedAttributes


    attributes.get
  }

  private def getAttributesWithDefaults(taskAttributes: Map[String, WdlValue], workflowOptions: Option[WorkflowOptions]): Map[String, WdlValue] = {
    val defaultAttributes = workflowOptions match {
      case Some(options) => defaultValues(options)
      case None => defaultValues
    }

    // Merge the two maps, prefer values from the latter
    defaultAttributes ++ taskAttributes
  }

  def validateKeys(keys: Set[String], backendType: BackendType): Try[Set[String]] = {
    /**
      *  Warn if keys are found that are unsupported on this backend.  This is not necessarily an error, at this point
      *  the keys are known to be supported on at least one backend other than this.
      */
    val unsupported = unsupportedKeys(keys, backendType)
    logMessageForUnsupportedKeys(unsupported, backendType) foreach log.warn

    // Return only supported keys
    validateRequiredKeys(keys, backendType) match {
      case scalaz.Success(x) => scala.util.Success(keys)
      case scalaz.Failure(e) => scala.util.Failure(new IllegalArgumentException() with ThrowableWithErrors {
        val message = "RuntimeAttributes are invalid."
        val errors = e
      })
    }
  }

  def logMessageForUnsupportedKeys(unsupportedKeys: Set[String], backendType: BackendType) = {
    if (unsupportedKeys.isEmpty) None
    else Option(s"Found unsupported keys for backend '$backendType': " + unsupportedKeys.toSeq.sorted.mkString(", "))
  }

  def unsupportedKeys(keys: Iterable[String], backendType: BackendType): Set[String] = {
    val supportedKeys = backendType.supportedKeys map { _.key }
    keys.toSet -- supportedKeys
  }

  val LocalDiskName = "local-disk"
  val LocalizationDisk = LocalDisk(LocalDiskName, DiskType.SSD, sizeGb=10)

  private case class ValidKeyType(key: String, validTypes: Set[WdlType])
  private val keys = Set(ValidKeyType("cpu", Set(WdlIntegerType)),
                        ValidKeyType("disks", Set(WdlStringType, WdlArrayType(WdlStringType))),
                        ValidKeyType("docker", Set(WdlStringType)),
                        ValidKeyType("zones", Set(WdlStringType, WdlArrayType(WdlStringType))),
                        ValidKeyType("continueOnReturnCode", Set(WdlBooleanType, WdlArrayType(WdlIntegerType))),
                        ValidKeyType("failOnStderr", Set(WdlBooleanType)),
                        ValidKeyType("preemptible", Set(WdlIntegerType)),
                        ValidKeyType("memory", Set(WdlStringType)))

  private val defaultValues = Map(
    "cpu" -> WdlInteger(1),
    "disks" -> WdlString(LocalizationDisk.toString),
    "zones" -> WdlString("us-central1-a"),
    "continueOnReturnCode" -> WdlInteger(0),
    "failOnStderr" -> WdlBoolean.False,
    "preemptible" -> WdlInteger(0),
    "memory" -> WdlString("2GB")
  )

  private def defaultValues(workflowOptions: WorkflowOptions): Map[String, WdlValue] = {
    def fromWorkflowOptions(k: ValidKeyType) = {
      val value = workflowOptions.getDefaultRuntimeOption(k.key).get
      val coercedValue = k.validTypes map { _.coerceRawValue(value) } find { _.isSuccess } map { _.get } getOrElse {
        log.warn(s"Could not coerce $value for runtime attribute ${k.key} to any of the supported target types: ${k.validTypes.mkString(",")}. Using default value instead.")
        defaultValues.get(k.key).get
      }
      coercedValue
    }

    keys.collect({
      case k if workflowOptions.getDefaultRuntimeOption(k.key).isSuccess => k.key -> fromWorkflowOptions(k)
      case k if defaultValues.contains(k.key) => k.key -> defaultValues.get(k.key).get
    }).toMap
  }

  val defaults: CromwellRuntimeAttributes = {
    // .get below because we assume that the default values are valid
    // There's a unit test in CromwellRuntimeAttributeSpec to test these defaults
    validateRuntimeAttributes(defaultValues).get
  }

  private def validateRuntimeAttributes(attributes: Map[String, WdlValue]): Try[CromwellRuntimeAttributes] = {
    val attributeMap = attributes.collect({ case (k, v) if Try(RuntimeKey.from(k)).isSuccess => RuntimeKey.from(k) -> v })
    val docker = validateDocker(attributeMap.get(DOCKER))
    val zones = validateZone(attributeMap.get(ZONES))
    val failOnStderr = validateFailOnStderr(attributeMap.get(FAIL_ON_STDERR))
    val continueOnReturnCode = validateContinueOnReturnCode(attributeMap.get(CONTINUE_ON_RETURN_CODE))
    val cpu = validateCpu(attributeMap.get(CPU))
    val disks = validateLocalDisks(attributeMap.get(DISKS))
    val preemptible = validatePreemptible(attributeMap.get(PREEMPTIBLE))
    val memory = RuntimeAttributes.validateMemoryValue(attributeMap.getOrElse(MEMORY, WdlString(s"${defaults.memoryGB} GB"))) match {
      case scala.util.Success(x) => x.to(MemoryUnit.GB).amount.successNel
      case scala.util.Failure(x) => x.getMessage.failureNel
    }

    (docker |@| zones |@| failOnStderr |@| continueOnReturnCode |@| cpu |@| preemptible |@| disks |@| memory) {
      new CromwellRuntimeAttributes(attributes, _, _, _, _, _, _, _, _)
    } match {
      case Success(x) => scala.util.Success(x)
      case Failure(nel) => scala.util.Failure(new IllegalArgumentException(nel.list.mkString("\n")))
    }
  }

  private def validateCpu(cpu: Option[WdlValue]): ErrorOr[Int] = {
    cpu.map(validateInt).getOrElse(defaults.cpu.toInt.successNel)
  }

  private def validatePreemptible(preemptible: Option[WdlValue]): ErrorOr[Int] = {
    preemptible.map(validateInt).getOrElse(defaults.preemptible.successNel)
  }

  private def validateDocker(docker: Option[WdlValue]): ErrorOr[Option[String]] = {
    docker match {
      case Some(WdlString(s)) => Some(s).successNel
      case None => None.successNel
      case _ => s"Expecting ${DOCKER.key} runtime attribute to a String".failureNel
    }
  }

  private def validateZone(zoneValue: Option[WdlValue]): ErrorOr[Vector[String]] = {
    zoneValue match {
      case Some(WdlString(s)) => s.split("\\s+").toVector.successNel
      case Some(WdlArray(wdlType, value)) if wdlType.memberType == WdlStringType =>
        value.map(_.valueString).toVector.successNel
      case Some(_) => s"Expecting ${ZONES.key} runtime attribute to be either a whitespace separated String or an Array[String]".failureNel
      case None => defaults.zones.toVector.successNel
    }
  }

  private def validateFailOnStderr(value: Option[WdlValue]): ErrorOr[Boolean] = {
    value match {
      case Some(WdlBoolean(b)) => b.successNel
      case Some(WdlString(s)) if s.toLowerCase == "true" => true.successNel
      case Some(WdlString(s)) if s.toLowerCase == "false" => false.successNel
      case Some(_) => s"Expecting ${FAIL_ON_STDERR.key} runtime attribute to be a Boolean or a String with values of 'true' or 'false'".failureNel
      case None => defaults.failOnStderr.successNel
    }
  }

  private def validateContinueOnReturnCode(value: Option[WdlValue]): ErrorOr[ContinueOnReturnCode] = {
    value match {
      case Some(b: WdlBoolean) => ContinueOnReturnCodeFlag(b.value).successNel
      case Some(WdlString(s)) if s.toLowerCase == "true" => ContinueOnReturnCodeFlag(true).successNel
      case Some(WdlString(s)) if s.toLowerCase == "false" => ContinueOnReturnCodeFlag(false).successNel
      case Some(WdlInteger(i)) => ContinueOnReturnCodeSet(Set(i)).successNel
      case Some(WdlArray(wdlType, seq)) =>
        val nels = seq map validateInt
        val defaultReturnCodeNel = Set.empty[Int].successNel[String]
        nels.foldLeft(defaultReturnCodeNel)((acc, v) => (acc |@| v) { (a, v) => a + v }) map ContinueOnReturnCodeSet
      case Some(_) => s"Expecting ${CONTINUE_ON_RETURN_CODE.key} runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]".failureNel
      case None => defaults.continueOnReturnCode.successNel
    }
  }

  private def validateLocalDisks(value: Option[WdlValue]): ErrorOr[Seq[Disk]] = {
    val nels = value match {
      case Some(WdlString(s)) => s.split(",\\s*").toSeq.map(validateLocalDisk)
      case Some(WdlArray(wdlType, seq)) if wdlType.memberType == WdlStringType =>
        seq.map(_.valueString).map(validateLocalDisk)
      case Some(_) =>
        Seq(s"Expecting ${DISKS.key} runtime attribute to be a comma separated String or Array[String]".failureNel[Disk])
      case None => defaults.disks.map(_.successNel)
    }

    val emptyDiskNel = Vector.empty[Disk].successNel[String]
    val disksNel = nels.foldLeft(emptyDiskNel)((acc, v) => (acc |@| v) { (a, v) => a :+ v })
    disksNel map {
      case disks if disks.exists(_.getName == LocalDiskName) => disks
      case disks => disks :+ LocalizationDisk.toDisk
    }
  }

  private def validateLocalDisk(disk: String): ErrorOr[Disk] = {
    LocalDisk.parse(disk) match {
      case scala.util.Success(localDisk) => localDisk.toDisk.successNel
      case scala.util.Failure(ex) => ex.getMessage.failureNel
    }
  }

  private def validateRequiredKeys(keySet: Set[String],
                                   backendType: BackendType): ErrorOr[Unit] = {
    val requiredKeys = for {
      key <- RuntimeKey.values().toSet
      if key.isMandatory(backendType)
    } yield key.key

    val missingRequiredKeys = requiredKeys -- keySet
    if (missingRequiredKeys.isEmpty) ().successNel
    else {
      val missingKeyString = missingRequiredKeys.toSeq.sorted.mkString(", ")
      s"Missing required keys in runtime configuration for backend '$backendType': $missingKeyString".failureNel
    }
  }

  private def validateBoolean(value: WdlValue): ErrorOr[Boolean] = {
    WdlBooleanType.coerceRawValue(value) match {
      case scala.util.Success(WdlBoolean(b)) => b.successNel
      case _ => s"Could not coerce $value into an boolean".failureNel
    }
  }

  private def validateInt(value: WdlValue): ErrorOr[Int] = {
    WdlIntegerType.coerceRawValue(value) match {
      case scala.util.Success(WdlInteger(i)) => i.intValue.successNel
      case _ => s"Could not coerce $value into an integer".failureNel
    }
  }

  implicit class EnhancedBackendType(val backendType: BackendType) extends AnyVal {
    def supportedKeys: Set[RuntimeKey] = for {
      key <- RuntimeKey.values().toSet
      if key.supports(backendType)
    } yield key
  }

  val DefaultReturnCodeNel = Set.empty[Int].successNel[String]
  val MemoryAndUnitPattern = """(\d+)\s*(\w+)""".r
}
