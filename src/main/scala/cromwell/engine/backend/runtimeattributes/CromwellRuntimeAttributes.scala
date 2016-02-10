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
import wdl4s.types.{WdlBooleanType, WdlIntegerType, WdlStringType}
import wdl4s.values._

import scala.language.postfixOps
import scala.util.Try
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

  def apply(wdlRuntimeAttributes: RuntimeAttributes, backendCall: BackendCall, workflowOptions: Option[WorkflowOptions]): CromwellRuntimeAttributes = {
    val attributes = for {
      keys <- validateKeys(wdlRuntimeAttributes, backendCall.backend.backendType)
      attributesFromTask <- TryUtil.sequenceMap(wdlRuntimeAttributes.evaluate(backendCall.lookupFunction(backendCall.locallyQualifiedInputs), backendCall.engineFunctions))
      attributesWithDefaults <- Try(getAttributesWithDefaults(attributesFromTask, workflowOptions))
      validatedAttributes <- validateRuntimeAttributes(attributesWithDefaults)
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

  def validateKeys(wdlRuntimeAttributes: RuntimeAttributes, backendType: BackendType): Try[RuntimeAttributes] = {
    /**
      *  Warn if keys are found that are unsupported on this backend.  This is not necessarily an error, at this point
      *  the keys are known to be supported on at least one backend other than this.
      */
    unsupportedKeys(wdlRuntimeAttributes.attrs.keys, backendType) foreach log.warn
    validateRequiredKeys(wdlRuntimeAttributes.attrs, backendType) match {
      case scalaz.Success(x) => scala.util.Success(wdlRuntimeAttributes)
      case scalaz.Failure(e) => scala.util.Failure(new IllegalArgumentException() with ThrowableWithErrors {
        val message = "RuntimeAttributes are invalid."
        val errors = e
      })
    }
  }

  def unsupportedKeys(keys: Iterable[String], backendType: BackendType): Option[String] = {
    val supportedKeys = backendType.supportedKeys map { _.key }
    val unsupportedKeys = keys.toSet -- supportedKeys

    if (unsupportedKeys.isEmpty) None
    else Option(s"Found unsupported keys for backend '$backendType': " + unsupportedKeys.toSeq.sorted.mkString(", "))
  }

  val LocalDiskName = "local-disk"
  val LocalizationDisk = LocalDisk(LocalDiskName, DiskType.SSD, sizeGb=10)

  private val keys = Set("cpu", "disks", "docker", "zones", "continueOnReturnCode", "failOnStderr", "preemptible", "memory")

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
    keys.collect({
      case k if workflowOptions.getDefaultRuntimeOption(k).isSuccess =>
        k -> WdlString(workflowOptions.getDefaultRuntimeOption(k).get)
      case k if defaultValues.contains(k) =>
        k -> defaultValues.get(k).get
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
      new CromwellRuntimeAttributes(_, _, _, _, _, _, _, _)
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

  private def validateRequiredKeys(attributeMap: Map[String, WdlExpression],
                                   backendType: BackendType): ErrorOr[Unit] = {
    val requiredKeys = for {
      key <- RuntimeKey.values().toSet
      if key.isMandatory(backendType)
    } yield key.key

    val missingRequiredKeys = requiredKeys -- attributeMap.keySet
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
