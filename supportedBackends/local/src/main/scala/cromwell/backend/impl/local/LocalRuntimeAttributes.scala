package cromwell.backend.impl.local

import cromwell.backend.validation.ContinueOnReturnCode
import cromwell.backend.validation.RuntimeAttributesDefault._
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation.RuntimeAttributesValidation._
import cromwell.core.WorkflowOptions
import lenthall.exception.MessageAggregation
import org.slf4j.Logger
import wdl4s.types._
import wdl4s.values.{WdlBoolean, WdlInteger, WdlValue}

import scalaz.Scalaz._
import scalaz._

object LocalRuntimeAttributes {
  val staticDefaults = Map(
    FailOnStderrKey -> WdlBoolean(false),
    ContinueOnReturnCodeKey -> WdlInteger(0)
  )

  val coercionMap: Map[String, Set[WdlType]] = Map (
    FailOnStderrKey -> Set[WdlType](WdlBooleanType),
    ContinueOnReturnCodeKey -> ContinueOnReturnCode.validWdlTypes,
    DockerKey -> Set(WdlStringType)
  )

  def apply(attrs: Map[String, WdlValue], options: WorkflowOptions, logger: Logger): LocalRuntimeAttributes = {
    // Fail now if some workflow options are specified but can't be parsed correctly
    val defaultFromOptions = workflowOptionsDefault(options, coercionMap).get
    val withDefaultValues = withDefaults(attrs, List(defaultFromOptions, staticDefaults))

    withDefaultValues.keySet.diff(coercionMap.keySet) map { k => s"Unrecognized runtime attribute key: $k" } foreach logger.warn

    val docker = validateDocker(withDefaultValues.get(DockerKey), None.successNel)
    val failOnStderr = validateFailOnStderr(withDefaultValues.get(FailOnStderrKey), noValueFoundFor(FailOnStderrKey))
    val continueOnReturnCode = validateContinueOnReturnCode(withDefaultValues.get(ContinueOnReturnCodeKey), noValueFoundFor(ContinueOnReturnCodeKey))
    (continueOnReturnCode |@| docker |@| failOnStderr) {
      new LocalRuntimeAttributes(_, _, _)
    } match {
      case Success(x) => x
      case Failure(nel) => throw new RuntimeException with MessageAggregation {
        override def exceptionContext: String = "Runtime attribute validation failed"
        override def errorMessages: Traversable[String] = nel.list
      }
    }
  }
}

case class LocalRuntimeAttributes(continueOnReturnCode: ContinueOnReturnCode, dockerImage: Option[String], failOnStderr: Boolean) {
  lazy val asMap = Map[String, Any](
    ContinueOnReturnCodeKey -> continueOnReturnCode,
    FailOnStderrKey -> failOnStderr.toString
  ) ++ dockerImage.map(img => Map(DockerKey -> img)).getOrElse(Map.empty)
}
