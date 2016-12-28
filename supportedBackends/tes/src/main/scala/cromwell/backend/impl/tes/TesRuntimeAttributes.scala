package cromwell.backend.impl.tes


import cats.data.Validated.{Invalid, Valid}
import cats.syntax.cartesian._
import cats.syntax.validated._
import cromwell.backend.MemorySize
import cromwell.backend.validation.ContinueOnReturnCode
import cromwell.backend.validation.RuntimeAttributesDefault._
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation.RuntimeAttributesValidation._
import cromwell.core.WorkflowOptions
import lenthall.validation.ErrorOr._
import lenthall.exception.MessageAggregation
import wdl4s.types.{WdlBooleanType, WdlIntegerType, WdlStringType, WdlType}
import wdl4s.values.{WdlBoolean, WdlInteger, WdlString, WdlValue}


object TesRuntimeAttributes {

  val DockerWorkingDirKey = "dockerWorkingDir"
  val DiskKey = "disk"

  val staticDefaults = Map(
    FailOnStderrKey         -> WdlBoolean(false),
    ContinueOnReturnCodeKey -> WdlInteger(0),
    CpuKey                  -> WdlInteger(1),
    MemoryKey               -> WdlString("1 GB"),
    DiskKey                 -> WdlString("1 GB")
  )

  private[tes] val coercionMap: Map[String, Set[WdlType]] = Map(
    FailOnStderrKey         -> Set[WdlType](WdlBooleanType),
    ContinueOnReturnCodeKey -> ContinueOnReturnCode.validWdlTypes,
    DockerKey               -> Set(WdlStringType),
    DockerWorkingDirKey     -> Set(WdlStringType),
    CpuKey                  -> Set(WdlIntegerType),
    MemoryKey               -> Set(WdlStringType),
    DiskKey                 -> Set(WdlStringType)
  )

  def apply(attrs: Map[String, WdlValue], options: WorkflowOptions): TesRuntimeAttributes = {
    // Fail now if some workflow options are specified but can't be parsed correctly
    val defaultFromOptions = workflowOptionsDefault(options, coercionMap).get
    val withDefaultValues = withDefaults(attrs, List(defaultFromOptions, staticDefaults))

    val docker               = validateDocker(withDefaultValues.get(DockerKey), noValueFoundFor(DockerKey))
    val dockerWorkingDir     = validateDockerWorkingDir(withDefaultValues.get(DockerWorkingDirKey), None.validNel)
    val failOnStderr         = validateFailOnStderr(withDefaultValues.get(FailOnStderrKey), noValueFoundFor(FailOnStderrKey))
    val continueOnReturnCode = validateContinueOnReturnCode(withDefaultValues.get(ContinueOnReturnCodeKey), noValueFoundFor(ContinueOnReturnCodeKey))
    val cpu                  = validateCpu(withDefaultValues.get(CpuKey), noValueFoundFor(CpuKey))
    val memory               = validateMemory(withDefaultValues.get(MemoryKey), noValueFoundFor(MemoryKey))
    val disk                 = validateDisk(withDefaultValues.get(DiskKey), noValueFoundFor(DiskKey))

    (continueOnReturnCode |@| docker |@| dockerWorkingDir |@| failOnStderr |@| cpu |@| memory |@| disk) map {
      new TesRuntimeAttributes(_, _, _, _, _, _, _)
    } match {
      case Valid(x: TesRuntimeAttributes) => x
      case Invalid(nel) => throw new RuntimeException with MessageAggregation {
        override def exceptionContext: String = "Runtime attribute validation failed"
        override def errorMessages: Traversable[String] = nel.toList
      }
    }
  }


  private def validateDockerWorkingDir(dockerWorkingDir: Option[WdlValue], onMissingKey: => ErrorOr[Option[String]]): ErrorOr[Option[String]] = {
    dockerWorkingDir match {
      case Some(WdlString(s)) => Some(s).validNel
      case None => onMissingKey
      case _ => s"Expecting $DockerWorkingDirKey runtime attribute to be a String".invalidNel
    }
  }

  private def validateDisk(value: Option[WdlValue], onMissingKey: => ErrorOr[MemorySize]): ErrorOr[MemorySize] = {
    val diskWrongFormatMsg = s"Expecting $DiskKey runtime attribute to be an Integer or String with format '8 GB'. Exception: %s"

    value match {
      case Some(i: WdlInteger) => parseMemoryInteger(i)
      case Some(s: WdlString) => parseMemoryString(s)
      case Some(_) => String.format(diskWrongFormatMsg, "Not supported WDL type value").invalidNel
      case None => onMissingKey
    }
  }
}

case class TesRuntimeAttributes(continueOnReturnCode: ContinueOnReturnCode,
                                dockerImage: Option[String],
                                dockerWorkingDir: Option[String],
                                failOnStderr: Boolean,
                                cpu: Int,
                                memory: MemorySize,
                                disk: MemorySize)
