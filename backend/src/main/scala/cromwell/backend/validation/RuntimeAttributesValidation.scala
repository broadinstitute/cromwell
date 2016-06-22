package cromwell.backend.validation

import cromwell.backend.MemorySize
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.core._
import org.slf4j.Logger
import wdl4s.parser.MemoryUnit
import wdl4s.types.WdlIntegerType
import wdl4s.values._

import scalaz.Scalaz._
import scalaz.{Failure, NonEmptyList}

object RuntimeAttributesValidation {
  val MemoryWrongAmountMsg = "Expecting %s runtime attribute value greater than 0 but got %s"
  val MemoryWrongFormatMsg = s"Expecting $MemoryKey runtime attribute to be an Integer or String with format '8 GB'. Exception: %s"

  def warnUnrecognized(actual: Set[String], expected: Set[String], logger: Logger) = {
    val unrecognized = actual.diff(expected).mkString(", ")
    if(unrecognized.nonEmpty) logger.warn(s"Unrecognized runtime attribute keys: $unrecognized")
  }

  def validateDocker(docker: Option[WdlValue], onMissingKey: => ErrorOr[Option[String]]): ErrorOr[Option[String]] = {
    docker match {
      case Some(WdlString(s)) => Some(s).successNel
      case None => onMissingKey
      case _ => s"Expecting $DockerKey runtime attribute to be a String".failureNel
    }
  }

  def validateFailOnStderr(value: Option[WdlValue], onMissingKey: => ErrorOr[Boolean]): ErrorOr[Boolean] = {
    value match {
      case Some(WdlBoolean(b)) => b.successNel
      case Some(WdlString(s)) if s.toLowerCase == "true" => true.successNel
      case Some(WdlString(s)) if s.toLowerCase == "false" => false.successNel
      case Some(_) => s"Expecting $FailOnStderrKey runtime attribute to be a Boolean or a String with values of 'true' or 'false'".failureNel
      case None => onMissingKey
    }
  }

  def validateContinueOnReturnCode(value: Option[WdlValue], onMissingKey: => ErrorOr[ContinueOnReturnCode]): ErrorOr[ContinueOnReturnCode] = {
    val failureWithMessage = s"Expecting $ContinueOnReturnCodeKey runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]".failureNel
    value match {
      case Some(b: WdlBoolean) => ContinueOnReturnCodeFlag(b.value).successNel
      case Some(WdlString(s)) if s.toLowerCase == "true" => ContinueOnReturnCodeFlag(true).successNel
      case Some(WdlString(s)) if s.toLowerCase == "false" => ContinueOnReturnCodeFlag(false).successNel
      case Some(WdlInteger(i)) => ContinueOnReturnCodeSet(Set(i)).successNel
      case Some(WdlArray(wdlType, seq)) =>
        val nels: Seq[ErrorOr[Int]] = seq map validateInt
        val nrFailures: Int = nels count { case p => p.isInstanceOf[Failure[NonEmptyList[String]]] }
        if (nrFailures == 0) {
          val defaultReturnCodeNel = Set.empty[Int].successNel[String]
          nels.foldLeft(defaultReturnCodeNel)((acc, v) => (acc |@| v) { (a, v) => a + v }) map ContinueOnReturnCodeSet
        }
        else failureWithMessage
      case Some(_) => failureWithMessage
      case None => onMissingKey
    }
  }

  def validateInt(value: WdlValue): ErrorOr[Int] = {
    WdlIntegerType.coerceRawValue(value) match {
      case scala.util.Success(WdlInteger(i)) => i.intValue.successNel
      case _ => s"Could not coerce $value into an integer".failureNel
    }
  }

  def validateMemory(value: Option[WdlValue], onMissingKey: => ErrorOr[MemorySize]): ErrorOr[MemorySize] = {
    value match {
      case Some(i: WdlInteger) => parseMemoryInteger(i)
      case Some(s: WdlString) => parseMemoryString(s)
      case Some(_) => String.format(MemoryWrongFormatMsg, "Not supported WDL type value").failureNel
      case None => onMissingKey
    }
  }

  def parseMemoryString(s: WdlString): ErrorOr[MemorySize] = {
    MemorySize.parse(s.valueString) match {
      case scala.util.Success(x: MemorySize) =>
        if (x.amount <= 0)
          String.format(MemoryWrongAmountMsg, MemoryKey, x.amount.toString()).failureNel
        else
          x.to(MemoryUnit.GB).successNel
      case scala.util.Failure(t) => String.format(MemoryWrongFormatMsg, t.getMessage).failureNel
    }
  }

  def parseMemoryInteger(i: WdlInteger): ErrorOr[MemorySize] = {
    if (i.value <= 0)
      String.format(MemoryWrongAmountMsg, MemoryKey, i.value.toString()).failureNel
    else
      MemorySize(i.value.toDouble, MemoryUnit.Bytes).to(MemoryUnit.GB).successNel
  }

  def validateCpu(cpu: Option[WdlValue], onMissingKey: => ErrorOr[Int]): ErrorOr[Int] = {
    val cpuValidation = cpu.map(validateInt).getOrElse(onMissingKey)
    cpuValidation match {
      case scalaz.Success(i) =>
        if (i <= 0)
          s"Expecting $CpuKey runtime attribute value greater than 0".failureNel
        else
          i.successNel
      case scalaz.Failure(f) =>
        s"Expecting $CpuKey runtime attribute to be an Integer".failureNel
    }
  }

}
