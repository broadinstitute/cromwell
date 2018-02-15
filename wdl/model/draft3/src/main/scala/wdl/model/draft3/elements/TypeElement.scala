package wdl.model.draft3.elements

import wom.types.WomPrimitiveType

import scala.language.implicitConversions

sealed trait TypeElement extends LanguageElement

case class PrimitiveTypeElement(primitiveType: WomPrimitiveType) extends TypeElement

case class ArrayTypeElement(inner: TypeElement) extends TypeElement
case class MapTypeElement(keyType: TypeElement, valueType: TypeElement ) extends TypeElement
case class OptionalTypeElement(maybeType: TypeElement) extends TypeElement
case class PairTypeElement(leftType: TypeElement, rightType: TypeElement) extends TypeElement
case class StructTypeElement(structName: String) extends TypeElement
case object ObjectTypeElement extends TypeElement
