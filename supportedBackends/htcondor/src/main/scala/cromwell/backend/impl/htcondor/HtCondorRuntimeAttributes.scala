package cromwell.backend.impl.htcondor

import cromwell.backend.MemorySize
import cromwell.backend.validation.ContinueOnReturnCode
import cromwell.backend.validation.RuntimeAttributesDefault._
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation.RuntimeAttributesValidation._
import cromwell.core._
import lenthall.exception.MessageAggregation
import wdl4s.types.{WdlIntegerType, WdlStringType, WdlBooleanType, WdlType}
import wdl4s.values.{WdlString, WdlBoolean, WdlInteger, WdlValue}

import scalaz.Scalaz._
import scalaz._

object HtCondorRuntimeAttributes {
  private val FailOnStderrDefaultValue = false
  private val ContinueOnRcDefaultValue = 0
  private val CpuDefaultValue = 1
  private val MemoryDefaultValue = "0.512 GB"
  private val DisksDefaultValue = "1.024 GB"

  val DockerWorkingDirKey = "dockerWorkingDir"
  val DockerOutputDirKey = "dockerOutputDir"
  val DiskKey = "disk"

  val staticDefaults = Map(
    FailOnStderrKey -> WdlBoolean(FailOnStderrDefaultValue),
    ContinueOnReturnCodeKey -> WdlInteger(ContinueOnRcDefaultValue),
    CpuKey -> WdlInteger(CpuDefaultValue),
    MemoryKey -> WdlString(MemoryDefaultValue),
    DiskKey -> WdlString(DisksDefaultValue)
  )

  val coercionMap: Map[String, Set[WdlType]] = Map (
    FailOnStderrKey -> Set[WdlType](WdlBooleanType),
    ContinueOnReturnCodeKey -> ContinueOnReturnCode.validWdlTypes,
    DockerKey -> Set(WdlStringType),
    DockerWorkingDirKey -> Set(WdlStringType),
    DockerOutputDirKey -> Set(WdlStringType),
    CpuKey -> Set(WdlIntegerType),
    MemoryKey -> Set(WdlStringType),
    DiskKey -> Set(WdlStringType)
  )

  def apply(attrs: Map[String, WdlValue], options: WorkflowOptions): HtCondorRuntimeAttributes = {
    // Fail now if some workflow options are specified but can't be parsed correctly
    val defaultFromOptions = workflowOptionsDefault(options, coercionMap).get
    val withDefaultValues = withDefaults(attrs, List(defaultFromOptions, staticDefaults))

    val docker = validateDocker(withDefaultValues.get(DockerKey), None.successNel)
    val dockerWorkingDir = validateDockerWorkingDir(withDefaultValues.get(DockerWorkingDirKey), None.successNel)
    val dockerOutputDir = validateDockerOutputDir(withDefaultValues.get(DockerOutputDirKey), None.successNel)
    val failOnStderr = validateFailOnStderr(withDefaultValues.get(FailOnStderrKey), noValueFoundFor(FailOnStderrKey))
    val continueOnReturnCode = validateContinueOnReturnCode(withDefaultValues.get(ContinueOnReturnCodeKey), noValueFoundFor(ContinueOnReturnCodeKey))
    val cpu = validateCpu(withDefaultValues.get(CpuKey), noValueFoundFor(CpuKey))
    val memory = validateMemory(withDefaultValues.get(MemoryKey), noValueFoundFor(MemoryKey))
    val disk = validateDisk(withDefaultValues.get(DiskKey), noValueFoundFor(DiskKey))

    (continueOnReturnCode |@| docker |@| dockerWorkingDir |@| dockerOutputDir |@| failOnStderr |@| cpu |@| memory |@| disk) {
      new HtCondorRuntimeAttributes(_, _, _, _, _, _, _, _)
    } match {
      case Success(x) => x
      case Failure(nel) => throw new RuntimeException with MessageAggregation {
        override def exceptionContext: String = "Runtime attribute validation failed"
        override def errorMessages: Traversable[String] = nel.list
      }
    }
  }

  private def validateDockerWorkingDir(dockerWorkingDir: Option[WdlValue], onMissingKey: => ErrorOr[Option[String]]): ErrorOr[Option[String]] = {
    dockerWorkingDir match {
      case Some(WdlString(s)) => Some(s).successNel
      case None => onMissingKey
      case _ => s"Expecting $DockerWorkingDirKey runtime attribute to be a String".failureNel
    }
  }

  private def validateDockerOutputDir(dockerOutputDir: Option[WdlValue], onMissingKey: => ErrorOr[Option[String]]): ErrorOr[Option[String]] = {
    dockerOutputDir match {
      case Some(WdlString(s)) => Some(s).successNel
      case None => onMissingKey
      case _ => s"Expecting $DockerOutputDirKey runtime attribute to be a String".failureNel
    }
  }

  private def validateDisk(value: Option[WdlValue], onMissingKey: => ErrorOr[MemorySize]): ErrorOr[MemorySize] = {
    val diskWrongFormatMsg = s"Expecting $DiskKey runtime attribute to be an Integer or String with format '8 GB'. Exception: %s"

    value match {
      case Some(i: WdlInteger) => parseMemoryInteger(i)
      case Some(s: WdlString) => parseMemoryString(s)
      case Some(_) => String.format(diskWrongFormatMsg, "Not supported WDL type value").failureNel
      case None => onMissingKey
    }
  }
}

case class HtCondorRuntimeAttributes(continueOnReturnCode: ContinueOnReturnCode,
                                     dockerImage: Option[String],
                                     dockerWorkingDir: Option[String],
                                     dockerOutputDir: Option[String],
                                     failOnStderr: Boolean,
                                     cpu: Int,
                                     memory: MemorySize,
                                     disk: MemorySize)
