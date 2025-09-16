package cromwell.backend.validation

import akka.stream.scaladsl.Flow.identityTraversalBuilder.attributes
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import liquibase.pro.packaged.is
import wom.RuntimeAttributesKeys
import wom.types.{WomArrayType, WomStringType, WomType}
import wom.values._

/**
  * Validates the "container" runtime attribute as a `String` or `Seq[String]`, returning it as `Seq[String]`.
  * If the user supplies a plain `String`, return it as a list of one element.
  *
  * `instance` returns a validation that errors when no attribute is specified.
  *
  * There is no default, however `optional` can be used return the validated value as an `Option`, wrapped in a `Some`,
  * if present, or `None` if not found.
  *
  * As of WDL 1.1 this attribute is preferred over the deprecated `docker` attribute. Previous to 1.1, it is not
 * supported, and is removed from runtime attrs during WDL parsing so as not to interfere with `docker`.
  */
object ContainerValidation {
  lazy val instance: RuntimeAttributesValidation[Seq[String]] = new ContainerValidation
  lazy val optional: OptionalRuntimeAttributesValidation[Seq[String]] = instance.optional
}

class ContainerValidation extends RuntimeAttributesValidation[Seq[String]] {
  override val key: String = RuntimeAttributesKeys.ContainerKey

  override def coercion: Iterable[WomType] = Set(WomStringType, WomArrayType(WomStringType))

  override def usedInCallCaching: Boolean = true

  override protected def missingValueMessage: String = "Can't find an attribute value for key container"

  override protected def invalidValueMessage(value: WomValue): String = super.missingValueMessage

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Seq[String]]] = {
    case WomString(value) => value.validNel.map(Seq(_))
    case WomArray(womType, values) if womType.memberType == WomStringType =>
      values.map(_.valueString).toSeq.validNel
  }

  // Before doing normal validation of this attribute, ensure that 'docker' is not also present. These two
  // different ways of specifying the container image are mutually exclusive.
  override def validate(values: Map[String, WomValue]): ErrorOr[Seq[String]] =
    if (values.contains(RuntimeAttributesKeys.DockerKey) && values.contains(RuntimeAttributesKeys.ContainerKey)) {
      s"Must provide only one of '${RuntimeAttributesKeys.DockerKey}' and '${RuntimeAttributesKeys.ContainerKey}' runtime attributes.".invalidNel
    } else super.validate(values)
}
