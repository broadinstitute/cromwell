package cromwell.engine.backend.runtimeattributes

import cromwell.backend.impl.jes.io.{JesAttachedDisk, JesWorkingDisk}
import cromwell.backend.validation.{ContinueOnReturnCode, RuntimeAttributesValidation}
import cromwell.core.{ErrorOr, WorkflowOptions}
import cromwell.engine.backend.runtimeattributes.OldStyleRuntimeKey._
import cromwell.engine.backend.{BackendType, OldStyleBackendCallJobDescriptor}
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

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class CromwellRuntimeAttributes(attributes: Map[String, WdlValue],
                                     docker: Option[String],
                                     zones: Seq[String],
                                     failOnStderr: Boolean,
                                     continueOnReturnCode: ContinueOnReturnCode,
                                     cpu: Long,
                                     preemptible: Int,
                                     disks: Seq[JesAttachedDisk],
                                     memoryGB: Double,
                                     bootDiskSizeGb: Int)

object CromwellRuntimeAttributes {
  private val log = LoggerFactory.getLogger("RuntimeAttributes")

  def apply(wdlRuntimeAttributes: RuntimeAttributes, jobDescriptor: OldStyleBackendCallJobDescriptor, workflowOptions: Option[WorkflowOptions]): CromwellRuntimeAttributes = {
    val supportedKeys = jobDescriptor.backend.backendType.supportedKeys map { _.key }

    val evaluatedAttr = wdlRuntimeAttributes.attrs mapValues {
      _.evaluate(jobDescriptor.lookupFunction(jobDescriptor.locallyQualifiedInputs), jobDescriptor.callEngineFunctions)
    }

    val attributes = for {
      attributesFromTask <- TryUtil.sequenceMap(evaluatedAttr)
      attributesWithDefaults <- Try(getAttributesWithDefaults(attributesFromTask, workflowOptions))
      _ <- validateKeys(attributesWithDefaults.keySet, jobDescriptor.backend.backendType)
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
      * Warn if keys are found that are unsupported on this backend.  This is not necessarily an error, at this point
      * the keys are known to be supported on at least one backend other than this.
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

  val DefaultJesWorkingDiskString = s"${JesWorkingDisk.Name} 10 SSD"
  val DefaultJesWorkingDisk = JesAttachedDisk.parse(DefaultJesWorkingDiskString).get

  private case class ValidKeyType(key: String, validTypes: Set[WdlType])

  private val keys = Set(
    ValidKeyType("cpu", Set(WdlIntegerType)),
    ValidKeyType("disks", Set(WdlStringType, WdlArrayType(WdlStringType))),
    ValidKeyType("docker", Set(WdlStringType)),
    ValidKeyType("zones", Set(WdlStringType, WdlArrayType(WdlStringType))),
    ValidKeyType("continueOnReturnCode", Set(WdlBooleanType, WdlArrayType(WdlIntegerType))),
    ValidKeyType("failOnStderr", Set(WdlBooleanType)),
    ValidKeyType("preemptible", Set(WdlIntegerType)),
    ValidKeyType("memory", Set(WdlStringType)),
    ValidKeyType("bootDiskSizeGb", Set(WdlStringType)))

  private val defaultValues = Map(
    "cpu" -> WdlInteger(1),
    "disks" -> WdlString(DefaultJesWorkingDiskString),
    "zones" -> WdlString("us-central1-a"),
    "continueOnReturnCode" -> WdlInteger(0),
    "failOnStderr" -> WdlBoolean.False,
    "preemptible" -> WdlInteger(0),
    "memory" -> WdlString("2GB"),
    "bootDiskSizeGb" -> WdlInteger(10)
  )

  private def defaultValues(workflowOptions: WorkflowOptions): Map[String, WdlValue] = {
    def fromWorkflowOptions(k: ValidKeyType): WdlValue = {
      val value = workflowOptions.getDefaultRuntimeOption(k.key).get
      val coercedValue = k.validTypes map { _.coerceRawValue(value) } find { _.isSuccess } map {
        _.get
      } getOrElse {
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
    val attributeMap = attributes.collect({ case (k, v) if Try(OldStyleRuntimeKey.from(k)).isSuccess => OldStyleRuntimeKey.from(k) -> v })
    val docker = RuntimeAttributesValidation.validateDocker(attributeMap.get(DOCKER), None.successNel)
    val zones = validateZone(attributeMap.get(ZONES))
    val failOnStderr = RuntimeAttributesValidation.validateFailOnStderr(attributeMap.get(FAIL_ON_STDERR), defaults.failOnStderr.successNel)
    val continueOnReturnCode = RuntimeAttributesValidation.validateContinueOnReturnCode(attributeMap.get(CONTINUE_ON_RETURN_CODE),
      defaults.continueOnReturnCode.successNel)
    val cpu = RuntimeAttributesValidation.validateCpu(attributeMap.get(CPU), defaults.cpu.toInt.successNel)
    val disks = validateLocalDisks(attributeMap.get(DISKS))
    val preemptible = validatePreemptible(attributeMap.get(PREEMPTIBLE))
    val memory = RuntimeAttributesValidation.validateMemory(attributeMap.get(MEMORY),
      RuntimeAttributesValidation.parseMemoryString(WdlString(s"${defaults.memoryGB} GB"))) match {
      case Success(x) => x.to(MemoryUnit.GB).amount.successNel
      case Failure(f) => Failure(f)
    }
    val bootDiskSize = validateBootDisk(attributeMap.get(BOOT_DISK))

    (docker |@| zones |@| failOnStderr |@| continueOnReturnCode |@| cpu |@| preemptible |@| disks |@| memory |@| bootDiskSize) {
      new CromwellRuntimeAttributes(attributes, _, _, _, _, _, _, _, _, _)
    } match {
      case Success(x) => scala.util.Success(x)
      case Failure(nel) => scala.util.Failure(new IllegalArgumentException(nel.list.mkString("\n")))
    }
  }

  private def validatePreemptible(preemptible: Option[WdlValue]): ErrorOr[Int] = {
    preemptible.map(validateInt).getOrElse(defaults.preemptible.successNel)
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

  private def validateLocalDisks(value: Option[WdlValue]): ErrorOr[Seq[JesAttachedDisk]] = {
    val nels = value match {
      case Some(WdlString(s)) => s.split(",\\s*").toSeq.map(validateLocalDisk)
      case Some(WdlArray(wdlType, seq)) if wdlType.memberType == WdlStringType =>
        seq.map(_.valueString).map(validateLocalDisk)
      case Some(_) =>
        Seq(s"Expecting ${DISKS.key} runtime attribute to be a comma separated String or Array[String]".failureNel[JesAttachedDisk])
      case None => defaults.disks.map(_.successNel)
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

  private def validateBootDisk(diskSize: Option[WdlValue]): ErrorOr[Int] = diskSize match {
    case Some(x) if WdlIntegerType.isCoerceableFrom(x.wdlType) =>
      WdlIntegerType.coerceRawValue(x) match {
        case scala.util.Success(x: WdlInteger) => x.value.intValue.successNel
        case scala.util.Success(coercionhasgoneawry) => s"Coercion was expected to create an Integer but instead got $coercionhasgoneawry".failureNel
        case scala.util.Failure(t) => t.getMessage.failureNel
      }
    case None => defaults.bootDiskSizeGb.successNel
  }

  private def validateRequiredKeys(keySet: Set[String],
                                   backendType: BackendType): ErrorOr[Unit] = {
    val requiredKeys = for {
      key <- OldStyleRuntimeKey.values().toSet
      if key.isMandatory(backendType)
    } yield key.key

    val missingRequiredKeys = requiredKeys -- keySet
    if (missingRequiredKeys.isEmpty) ().successNel
    else {
      val missingKeyString = missingRequiredKeys.toSeq.sorted.mkString(", ")
      s"Missing required keys in runtime configuration for backend '$backendType': $missingKeyString".failureNel
    }
  }

  private def validateInt(value: WdlValue): ErrorOr[Int] = {
    WdlIntegerType.coerceRawValue(value) match {
      case scala.util.Success(WdlInteger(i)) => i.intValue.successNel
      case _ => s"Could not coerce $value into an integer".failureNel
    }
  }

  implicit class EnhancedBackendType(val backendType: BackendType) extends AnyVal {
    def supportedKeys: Set[OldStyleRuntimeKey] = for {
      key <- OldStyleRuntimeKey.values().toSet
      if key.supports(backendType)
    } yield key
  }
}
