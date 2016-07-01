package cromwell.backend.impl.spark

import cromwell.backend.validation.ContinueOnReturnCode
import cromwell.backend.validation.RuntimeAttributesDefault._
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation.RuntimeAttributesValidation._
import cromwell.core._
import lenthall.exception.MessageAggregation
import wdl4s.types.{WdlBooleanType, WdlStringType, WdlType}
import wdl4s.values.{WdlBoolean, WdlInteger, WdlString, WdlValue}

import scalaz.Scalaz._
import scalaz._

object SparkRuntimeAttributes {
  private val FailOnStderrDefaultValue = false
  private val ContinueOnRcDefaultValue = 0
  private val DockerWorkingDirKey = "dockerWorkingDir"
  private val DockerOutputDirKey = "dockerOutputDir"

  val staticDefaults = Map(
    FailOnStderrKey -> WdlBoolean(FailOnStderrDefaultValue),
    ContinueOnReturnCodeKey -> WdlInteger(ContinueOnRcDefaultValue)
  )

  val coercionMap: Map[String, Set[WdlType]] = Map(
    FailOnStderrKey -> Set[WdlType](WdlBooleanType),
    ContinueOnReturnCodeKey -> ContinueOnReturnCode.validWdlTypes,
    DockerKey -> Set(WdlStringType),
    DockerWorkingDirKey -> Set(WdlStringType),
    DockerOutputDirKey -> Set(WdlStringType)
  )

  def apply(attrs: Map[String, WdlValue], options: WorkflowOptions): SparkRuntimeAttributes = {
    // Fail now if some workflow options are specified but can't be parsed correctly
    val defaultFromOptions = workflowOptionsDefault(options, coercionMap).get
    val withDefaultValues = withDefaults(attrs, List(defaultFromOptions, staticDefaults))

    val docker = validateDocker(withDefaultValues.get(DockerKey), None.successNel)
    val dockerWorkingDir = validateDockerWorkingDir(withDefaultValues.get(DockerWorkingDirKey), None.successNel)
    val dockerOutputDir = validateDockerOutputDir(withDefaultValues.get(DockerOutputDirKey), None.successNel)
    val failOnStderr = validateFailOnStderr(withDefaultValues.get(FailOnStderrKey), noValueFoundFor(FailOnStderrKey))
    val continueOnReturnCode = validateContinueOnReturnCode(withDefaultValues.get(ContinueOnReturnCodeKey), noValueFoundFor(ContinueOnReturnCodeKey))


    (continueOnReturnCode |@| docker |@| dockerWorkingDir |@| dockerOutputDir |@| failOnStderr) {
      new SparkRuntimeAttributes(_, _, _, _, _)
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
}


case class SparkRuntimeAttributes(continueOnReturnCode: ContinueOnReturnCode,
                                  dockerImage: Option[String],
                                  dockerWorkingDir: Option[String],
                                  dockerOutputDir: Option[String],
                                  failOnStderr: Boolean)

