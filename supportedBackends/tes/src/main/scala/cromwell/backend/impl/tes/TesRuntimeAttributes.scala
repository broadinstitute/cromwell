package cromwell.backend.impl.tes


import cats.data.Validated.{Invalid, Valid}
import cats.syntax.cartesian._
import cats.syntax.validated._
import cromwell.backend.{BackendJobDescriptor, MemorySize}
import cromwell.backend.validation.ContinueOnReturnCode
import cromwell.backend.validation.RuntimeAttributesDefault._
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation.RuntimeAttributesValidation._
import cromwell.core._
import cromwell.core.ErrorOr._
import lenthall.exception.MessageAggregation
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.types.{WdlBooleanType, WdlIntegerType, WdlStringType, WdlType}
import wdl4s.util.TryUtil
import wdl4s.values.{WdlBoolean, WdlInteger, WdlString, WdlValue}


object TesRuntimeAttributes {
  private val FailOnStderrDefaultValue = false
  private val ContinueOnRcDefaultValue = 0
  private val CpuDefaultValue = 1
  private val MemoryDefaultValue = "1 GB"
  private val DisksDefaultValue = "1 GB"

  val DockerWorkingDirKey = "dockerWorkingDir"
  val DiskKey = "disk"

  val staticDefaults = Map(
    FailOnStderrKey         -> WdlBoolean(FailOnStderrDefaultValue),
    ContinueOnReturnCodeKey -> WdlInteger(ContinueOnRcDefaultValue),
    CpuKey                  -> WdlInteger(CpuDefaultValue),
    MemoryKey               -> WdlString(MemoryDefaultValue),
    DiskKey                 -> WdlString(DisksDefaultValue)
  )

  private[tes] val coercionMap: Map[String, Set[WdlType]] = Map (
    FailOnStderrKey         -> Set[WdlType](WdlBooleanType),
    ContinueOnReturnCodeKey -> ContinueOnReturnCode.validWdlTypes,
    DockerKey               -> Set(WdlStringType),
    DockerWorkingDirKey     -> Set(WdlStringType),
    CpuKey                  -> Set(WdlIntegerType),
    MemoryKey               -> Set(WdlStringType),
    DiskKey                 -> Set(WdlStringType)
  )

  def apply(attrs: Map[String, WdlValue]): TesRuntimeAttributes = {

    val docker               = validateDocker(attrs.get(DockerKey), noValueFoundFor(DockerKey))
    val dockerWorkingDir     = validateDockerWorkingDir(attrs.get(DockerWorkingDirKey), None.validNel)
    val failOnStderr         = validateFailOnStderr(attrs.get(FailOnStderrKey), noValueFoundFor(FailOnStderrKey))
    val continueOnReturnCode = validateContinueOnReturnCode(attrs.get(ContinueOnReturnCodeKey), noValueFoundFor(ContinueOnReturnCodeKey))
    val cpu                  = validateCpu(attrs.get(CpuKey), noValueFoundFor(CpuKey))
    val memory               = validateMemory(attrs.get(MemoryKey), noValueFoundFor(MemoryKey))
    val disk                 = validateDisk(attrs.get(DiskKey), noValueFoundFor(DiskKey))

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

  def fromJobDescriptor(jobDescriptor: BackendJobDescriptor, expressionFunctions: WdlStandardLibraryFunctions): TesRuntimeAttributes = {
    val lookup = jobDescriptor.inputs.apply _
    val evaluateAttrs = jobDescriptor.call.task
      .runtimeAttributes
      .attrs
      .mapValues(_.evaluate(lookup, expressionFunctions))
    // Fail the call if runtime attributes can't be evaluated
    val runtimeMap = TryUtil.sequenceMap(evaluateAttrs, "Runtime attributes evaluation").get
    // Fail now if some workflow options are specified but can't be parsed correctly
    val options = jobDescriptor.workflowDescriptor.workflowOptions
    val defaultFromOptions = workflowOptionsDefault(options, TesRuntimeAttributes.coercionMap).get
    val withDefaultValues = withDefaults(runtimeMap, List(defaultFromOptions, TesRuntimeAttributes.staticDefaults))
    TesRuntimeAttributes(withDefaultValues)
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
