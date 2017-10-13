package cromwell.backend.impl.spark

import cats.data.Validated.{Invalid, Valid}
import cats.syntax.apply._
import cats.syntax.validated._
import cromwell.backend.MemorySize
import cromwell.backend.validation.RuntimeAttributesDefault._
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation.RuntimeAttributesValidation._
import cromwell.core._
import lenthall.exception.MessageAggregation
import lenthall.validation.ErrorOr._
import wom.types._
import wom.values._

object SparkRuntimeAttributes {
  private val FailOnStderrDefaultValue = false
  private val ExecutorCoresDefaultValue = 1
  private val ExecutorMemoryDefaultValue = "1 GB"

  val ExecutorCoresKey = "executorCores"
  val ExecutorMemoryKey = "executorMemory"
  val AppMainClassKey = "appMainClass"
  //Specific to cluster mode
  val NumberOfExecutorsKey = "numberOfExecutors"
  val AdditionalArgsKey = "additionalArgs"

  val staticDefaults = Map(
    FailOnStderrKey -> WdlBoolean(FailOnStderrDefaultValue),
    ExecutorCoresKey -> WdlInteger(ExecutorCoresDefaultValue),
    ExecutorMemoryKey -> WdlString(ExecutorMemoryDefaultValue)
  )

  val coercionMap: Map[String, Set[WdlType]] = Map(
    FailOnStderrKey -> Set[WdlType](WdlBooleanType),
    ExecutorCoresKey -> Set(WdlIntegerType),
    ExecutorMemoryKey -> Set(WdlStringType),
    AppMainClassKey -> Set(WdlStringType),
    NumberOfExecutorsKey -> Set(WdlIntegerType),
    AdditionalArgsKey -> Set(WdlStringType)
  )

  def apply(attrs: Map[String, WdlValue], options: WorkflowOptions): SparkRuntimeAttributes = {
    // Fail now if some workflow options are specified but can't be parsed correctly
    val defaultFromOptions = workflowOptionsDefault(options, coercionMap).get
    val withDefaultValues = withDefaults(attrs, List(defaultFromOptions, staticDefaults))

    val failOnStderr = validateFailOnStderr(withDefaultValues.get(FailOnStderrKey), noValueFoundFor(FailOnStderrKey))

    val executorCores = validateCpu(withDefaultValues.get(ExecutorCoresKey), noValueFoundFor(ExecutorCoresKey))
    val executorMemory = validateMemory(withDefaultValues.get(ExecutorMemoryKey), noValueFoundFor(ExecutorMemoryKey))
    val numberOfExecutors = validateNumberOfExecutors(withDefaultValues.get(NumberOfExecutorsKey), None.validNel)
    val appMainClass = validateAppEntryPoint(withDefaultValues.get(AppMainClassKey), None.validNel)
    val additionalArgs = validateAdditionalArgs(withDefaultValues.get(AdditionalArgsKey), None.validNel)

    (executorCores, executorMemory, numberOfExecutors, appMainClass, additionalArgs, failOnStderr) mapN  {SparkRuntimeAttributes.apply} match {
      case Valid(x) => x
      case Invalid(nel) => throw new RuntimeException with MessageAggregation {
        override def exceptionContext: String = "Runtime attribute validation failed"
        override def errorMessages: Traversable[String] = nel.toList
      }
    }
  }

  private def validateNumberOfExecutors(numOfExecutors: Option[WdlValue], onMissingKey: => ErrorOr[Option[Int]]): ErrorOr[Option[Int]] = {
    numOfExecutors match {
      case Some(i: WdlInteger) => Option(i.value.intValue()).validNel
      case None => onMissingKey
      case _ => s"Expecting $NumberOfExecutorsKey runtime attribute to be an Integer".invalidNel
    }
  }

  private def validateAppEntryPoint(mainClass:  Option[WdlValue], onMissingKey: => ErrorOr[Option[String]]): ErrorOr[Option[String]] = {
    mainClass match {
      case Some(WdlString(s)) => Option(s).validNel
      case None => onMissingKey
      case _ => s"Expecting $AppMainClassKey runtime attribute to be a String".invalidNel
    }
  }

  private def validateAdditionalArgs(additionalArgs:  Option[WdlValue], onMissingKey: => ErrorOr[Option[String]]): ErrorOr[Option[String]] = {
    additionalArgs match {
      case Some(WdlString(s)) => Option(s).validNel
      case None => onMissingKey
      case _ => s"Expecting $AdditionalArgsKey runtime attribute to be a String".invalidNel
    }
  }
}

case class SparkRuntimeAttributes(executorCores: Int, executorMemory: MemorySize, numberOfExecutors: Option[Int],
                                  appMainClass: Option[String], additionalArgs: Option[String], failOnStderr: Boolean)
