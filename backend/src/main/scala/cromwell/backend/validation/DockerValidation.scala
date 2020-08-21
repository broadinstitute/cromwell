package cromwell.backend.validation

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wom.RuntimeAttributesKeys
import wom.values._

/**
  * Validates the "docker" runtime attribute as a String, returning it as `String`.
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * There is no default, however `optional` can be used return the validated value as an `Option`, wrapped in a `Some`,
  * if present, or `None` if not found.
  */
object DockerValidation {
  lazy val instance: RuntimeAttributesValidation[String] = new DockerValidation
  lazy val optional: OptionalRuntimeAttributesValidation[String] = instance.optional
}

class DockerValidation extends StringRuntimeAttributesValidation(RuntimeAttributesKeys.DockerKey) {
  override def usedInCallCaching: Boolean = true

  override protected def missingValueMessage: String = "Can't find an attribute value for key docker"

  override protected def invalidValueMessage(value: WomValue): String = super.missingValueMessage

  // NOTE: Docker's current test specs don't like WdlInteger, etc. auto converted to WdlString.
  override protected def validateValue: PartialFunction[WomValue, ErrorOr[String]] = {
    case WomString(value) => value.validNel
  }
}
