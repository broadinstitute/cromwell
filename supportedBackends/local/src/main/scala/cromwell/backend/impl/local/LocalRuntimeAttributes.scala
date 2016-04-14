package cromwell.backend.impl.local

import cromwell.backend.validation.{ContinueOnReturnCodeSet, ContinueOnReturnCode}
import cromwell.backend.validation.RuntimeAttributesKeys._
import wdl4s.values.WdlValue
import cromwell.backend.validation.RuntimeAttributesValidation._
import scalaz._
import Scalaz._
import lenthall.exception.MessageAggregation

object LocalRuntimeAttributes {
  val FailOnStderrDefaultValue = false
  val ContinueOnRcDefaultValue = 0

  def apply(attrs: Map[String, WdlValue]): LocalRuntimeAttributes = {

    // FIXME Duplicated from validation - need a way to share this
    val docker = validateDocker(attrs.get(Docker), () => None.successNel)
    val failOnStderr = validateFailOnStderr(attrs.get(FailOnStderr), () => FailOnStderrDefaultValue.successNel)
    val continueOnReturnCode = validateContinueOnReturnCode(attrs.get(ContinueOnReturnCode),
      () => ContinueOnReturnCodeSet(Set(ContinueOnRcDefaultValue)).successNel)
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

case class LocalRuntimeAttributes(continueOnReturnCode: ContinueOnReturnCode, dockerImage: Option[String], failOnStderr: Boolean)
