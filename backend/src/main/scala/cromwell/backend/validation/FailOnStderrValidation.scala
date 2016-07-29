package cromwell.backend.validation

import wdl4s.types.{WdlBooleanType, WdlStringType}
import wdl4s.values.{WdlBoolean, WdlString}

import scala.util.Try
import scalaz.Scalaz._

/**
  * Validates the "failOnStderr" runtime attribute as a Boolean or a String 'true' or 'false', returning the value as a
  * `Boolean`.
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * The default returns `false` when no attribute is specified.
  */
object FailOnStderrValidation {
  val key = RuntimeAttributesKeys.FailOnStderrKey

  lazy val instance = new FailOnStderrValidation

  lazy val default = instance.withDefault(WdlBoolean(false))

  lazy val optional = default.optional

  private[validation] val missingMessage =
    s"Expecting $key runtime attribute to be a Boolean or a String with values of 'true' or 'false'"
}

class FailOnStderrValidation extends RuntimeAttributesValidation[Boolean] {

  import FailOnStderrValidation._

  override def key = RuntimeAttributesKeys.FailOnStderrKey

  override def coercion = Seq(WdlBooleanType, WdlStringType)

  override protected def validateValue = {
    case WdlBoolean(value) => value.successNel
    case WdlString(value) if Try(value.toBoolean).isSuccess => value.toBoolean.successNel
  }

  override def validateExpression = {
    case _: WdlBoolean => true
    case WdlString(value) if Try(value.toBoolean).isSuccess => true
  }

  override protected def failureMessage = missingMessage
}
