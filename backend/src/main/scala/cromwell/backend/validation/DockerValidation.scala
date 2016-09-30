package cromwell.backend.validation

import cats.syntax.validated._
import wdl4s.types.WdlStringType
import wdl4s.values.WdlString

/**
  * Validates the "docker" runtime attribute as a String, returning it as `String`.
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * There is no default, however `optional` can be used return the validated value as an `Option`, wrapped in a `Some`,
  * if present, or `None` if not found.
  */
object DockerValidation {
  val key = RuntimeAttributesKeys.DockerKey

  lazy val instance = new DockerValidation

  lazy val optional = instance.optional

  private[validation] val missingMessage = s"Expecting $key runtime attribute to be a String"
}

class DockerValidation extends RuntimeAttributesValidation[String] {

  import DockerValidation._

  override def key = RuntimeAttributesKeys.DockerKey

  override def coercion = Seq(WdlStringType)

  override protected def validateValue = {
    case WdlString(value) => value.validNel
  }

  override protected def failureMessage = missingMessage
}
