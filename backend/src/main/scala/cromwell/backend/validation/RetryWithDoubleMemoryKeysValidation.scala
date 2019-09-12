package cromwell.backend.validation

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wom.RuntimeAttributesKeys
import wom.types.{WomArrayType, WomStringType, WomType}
import wom.values.{WomArray, WomValue}


object RetryWithDoubleMemoryKeysValidation {

  def instance(attributeName: String = RuntimeAttributesKeys.RetryWithDoubleMemoryKeys): RuntimeAttributesValidation[Vector[String]] =
    new RetryWithDoubleMemoryKeysValidation(attributeName)

  def optional(attributeName: String = RuntimeAttributesKeys.RetryWithDoubleMemoryKeys): OptionalRuntimeAttributesValidation[Vector[String]] =
    instance(attributeName).optional
}


class RetryWithDoubleMemoryKeysValidation(attributeName: String = RuntimeAttributesKeys.RetryWithDoubleMemoryKeys) extends RuntimeAttributesValidation[Vector[String]] {

  override def key: String = RuntimeAttributesKeys.RetryWithDoubleMemoryKeys

  override def coercion: Traversable[WomType] = Set(WomArrayType(WomStringType))

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Vector[String]]] = {
    case WomArray(womType, value) if womType.memberType == WomStringType => value.map(_.valueString).toVector.validNel
  }
}
