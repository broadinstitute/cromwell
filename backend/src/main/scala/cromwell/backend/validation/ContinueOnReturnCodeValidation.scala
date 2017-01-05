package cromwell.backend.validation

import cats.data.Validated.{Invalid, Valid}
import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import cromwell.backend.validation.RuntimeAttributesValidation._
import lenthall.validation.ErrorOr._
import wdl4s.types.{WdlArrayType, WdlIntegerType, WdlStringType}
import wdl4s.values.{WdlArray, WdlBoolean, WdlInteger, WdlString}

import scala.util.Try

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
    case WdlBoolean(value) => ContinueOnReturnCodeFlag(value).validNel
    case WdlString(value) if Try(value.toBoolean).isSuccess => ContinueOnReturnCodeFlag(value.toBoolean).validNel
    case WdlInteger(value) => ContinueOnReturnCodeSet(Set(value)).validNel
    case WdlArray(wdlType, seq) =>
      val errorOrInts: ErrorOr[List[Int]] = (seq.toList map validateInt).sequence[ErrorOr, Int]
      errorOrInts match {
        case Valid(ints) => ContinueOnReturnCodeSet(ints.toSet).validNel
        case Invalid(_) => failureWithMessage
      }
  }

  override def validateExpression = {
    case _: WdlBoolean => true
    case WdlString(value) if Try(value.toBoolean).isSuccess => true
    case _: WdlInteger => true
    case WdlArray(WdlArrayType(WdlStringType), elements) => elements.forall(validateInt(_).isValid)
    case WdlArray(WdlArrayType(WdlIntegerType), elements) => elements.forall(validateInt(_).isValid)
  }

  override protected def failureMessage = missingMessage
}
