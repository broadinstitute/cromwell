package cromwell.backend.impl.htcondor

import cromwell.backend.validation.ContinueOnReturnCode
import cromwell.backend.validation.RuntimeAttributesDefault._
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation.RuntimeAttributesValidation._
import cromwell.core.WorkflowOptions
import lenthall.exception.MessageAggregation
import wdl4s.types.{WdlStringType, WdlBooleanType, WdlType}
import wdl4s.values.{WdlBoolean, WdlInteger, WdlValue}

import scalaz.Scalaz._
import scalaz._

object HtCondorRuntimeAttributes {
  val FailOnStderrDefaultValue = false
  val ContinueOnRcDefaultValue = 0

  val staticDefaults = Map(
    FailOnStderrKey -> WdlBoolean(FailOnStderrDefaultValue),
    ContinueOnReturnCodeKey -> WdlInteger(ContinueOnRcDefaultValue)
  )

  val coercionMap: Map[String, Set[WdlType]] = Map (
    FailOnStderrKey -> Set[WdlType](WdlBooleanType),
    ContinueOnReturnCodeKey -> ContinueOnReturnCode.validWdlTypes,
    DockerKey -> Set(WdlStringType)
  )

  def apply(attrs: Map[String, WdlValue], options: WorkflowOptions): HtCondorRuntimeAttributes = {
    // Fail now if some workflow options are specified but can't be parsed correctly
    val defaultFromOptions = workflowOptionsDefault(options, coercionMap).get
    val withDefaultValues = withDefaults(attrs, List(defaultFromOptions, staticDefaults))

    val docker = validateDocker(withDefaultValues.get(DockerKey), None.successNel)
    val failOnStderr = validateFailOnStderr(withDefaultValues.get(FailOnStderrKey), noValueFoundFor(FailOnStderrKey))
    val continueOnReturnCode = validateContinueOnReturnCode(withDefaultValues.get(ContinueOnReturnCodeKey), noValueFoundFor(ContinueOnReturnCodeKey))
    (continueOnReturnCode |@| docker |@| failOnStderr) {
      new HtCondorRuntimeAttributes(_, _, _)
    } match {
      case Success(x) => x
      case Failure(nel) => throw new RuntimeException with MessageAggregation {
        override def exceptionContext: String = "Runtime attribute validation failed"
        override def errorMessages: Traversable[String] = nel.list
      }
    }
  }
}

case class HtCondorRuntimeAttributes(continueOnReturnCode: ContinueOnReturnCode, dockerImage: Option[String], failOnStderr: Boolean)
