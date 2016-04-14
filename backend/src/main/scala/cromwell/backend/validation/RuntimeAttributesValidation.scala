package cromwell.backend.validation

import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.core._
import wdl4s.types.WdlIntegerType
import wdl4s.values._

import scalaz.{NonEmptyList, Failure}
import scalaz.Scalaz._

object RuntimeAttributesValidation {

  def validateDocker(docker: Option[WdlValue], onMissingKey: () => ErrorOr[Option[String]]): ErrorOr[Option[String]] = {
    docker match {
      case Some(WdlString(s)) => Some(s).successNel
      case None => onMissingKey.apply()
      case _ => s"Expecting $Docker runtime attribute to be a String".failureNel
    }
  }

  def validateFailOnStderr(value: Option[WdlValue], onMissingKey: () => ErrorOr[Boolean]): ErrorOr[Boolean] = {
    value match {
      case Some(WdlBoolean(b)) => b.successNel
      case Some(WdlString(s)) if s.toLowerCase == "true" => true.successNel
      case Some(WdlString(s)) if s.toLowerCase == "false" => false.successNel
      case Some(_) => s"Expecting $FailOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'".failureNel
      case None => onMissingKey.apply()
    }
  }

  def validateContinueOnReturnCode(value: Option[WdlValue], onMissingKey: () => ErrorOr[ContinueOnReturnCode]): ErrorOr[ContinueOnReturnCode] = {
    val failureWithMessage = s"Expecting $ContinueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]".failureNel
    value match {
      case Some(b: WdlBoolean) => ContinueOnReturnCodeFlag(b.value).successNel
      case Some(WdlString(s)) if s.toLowerCase == "true" => ContinueOnReturnCodeFlag(true).successNel
      case Some(WdlString(s)) if s.toLowerCase == "false" => ContinueOnReturnCodeFlag(false).successNel
      case Some(WdlInteger(i)) => ContinueOnReturnCodeSet(Set(i)).successNel
      case Some(WdlArray(wdlType, seq)) =>
        val nels: Seq[ErrorOr[Int]] = seq map validateInt
        val nrFailures: Int = nels count { case p => p.isInstanceOf[Failure[NonEmptyList[String]]] }
        if(nrFailures == 0) {
          val defaultReturnCodeNel = Set.empty[Int].successNel[String]
          nels.foldLeft(defaultReturnCodeNel)((acc, v) => (acc |@| v) { (a, v) => a + v }) map ContinueOnReturnCodeSet
        }
        else {
          failureWithMessage
        }
      case Some(_) => failureWithMessage
      case None => onMissingKey.apply()
    }
  }

  def validateInt(value: WdlValue): ErrorOr[Int] = {
    WdlIntegerType.coerceRawValue(value) match {
      case scala.util.Success(WdlInteger(i)) => i.intValue.successNel
      case _ => s"Could not coerce $value into an integer".failureNel
    }
  }

}
