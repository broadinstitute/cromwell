package cromwell.backend.impl.htcondor

import cromwell.backend.validation.{ContinueOnReturnCodeSet, ContinueOnReturnCode}
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation.RuntimeAttributesValidation._
import lenthall.exception.MessageAggregation
import wdl4s.values.WdlValue

import scalaz._
import Scalaz._

object HtCondorRuntimeAttributes {
  val FailOnStderrDefaultValue = false
  val ContinueOnRcDefaultValue = 0

  def apply(attrs: Map[String, WdlValue]): HtCondorRuntimeAttributes = {
    val docker = validateDocker(attrs.get(Docker), None.successNel)
    val failOnStderr = validateFailOnStderr(attrs.get(FailOnStderr), FailOnStderrDefaultValue.successNel)
    val continueOnReturnCode = validateContinueOnReturnCode(attrs.get(ContinueOnReturnCode),
      ContinueOnReturnCodeSet(Set(ContinueOnRcDefaultValue)).successNel)
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
