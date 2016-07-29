package cromwell.backend.validation

import cromwell.backend.validation.RuntimeAttributesValidation._
import cromwell.core._
import wdl4s.types.{WdlArrayType, WdlStringType}
import wdl4s.values.{WdlArray, WdlBoolean, WdlInteger, WdlString}

import scala.util.Try
import scalaz.Scalaz._

/**
  * Validates the "continueOnReturnCode" runtime attribute a Boolean, a String 'true' or 'false', or an Array[Int],
  * returning the value as a instance of a class extending trait `ContinueOnReturnCode`.
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * The default returns a `ContinueOnReturnCodeSet(0)` when no attribute is specified.
  *
  * `optional` can be used return the validated value as an `Option`, wrapped in a `Some`, if present, or `None` if not
  * found.
  */
object ContinueOnReturnCodeValidation {
  val key = RuntimeAttributesKeys.ContinueOnReturnCodeKey

  lazy val instance = new ContinueOnReturnCodeValidation

  lazy val default = instance.withDefault(WdlInteger(0))

  lazy val optional = default.optional

  private[validation] val missingMessage =
    s"Expecting $key runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]"
}

class ContinueOnReturnCodeValidation extends RuntimeAttributesValidation[ContinueOnReturnCode] {

  import ContinueOnReturnCodeValidation._

  override def key = RuntimeAttributesKeys.ContinueOnReturnCodeKey

  override def coercion = ContinueOnReturnCode.validWdlTypes

  override def validateValue = {
    case WdlBoolean(value) => ContinueOnReturnCodeFlag(value).successNel
    case WdlString(value) if Try(value.toBoolean).isSuccess => ContinueOnReturnCodeFlag(value.toBoolean).successNel
    case WdlInteger(value) => ContinueOnReturnCodeSet(Set(value)).successNel
    case WdlArray(wdlType, seq) =>
      val errorOrInts: ErrorOr[List[Int]] = (seq.toList map validateInt).sequence[ErrorOr, Int]
      errorOrInts match {
        case scalaz.Success(ints) => ContinueOnReturnCodeSet(ints.toSet).successNel
        case scalaz.Failure(_) => failureWithMessage
      }
  }

  override def validateExpression = {
    case _: WdlBoolean => true
    case WdlString(value) if Try(value.toBoolean).isSuccess => true
    case _: WdlInteger => true
    case WdlArray(WdlArrayType(WdlStringType), elements) => elements.forall(validateInt(_).isSuccess)
  }

  override protected def failureMessage = missingMessage
}
