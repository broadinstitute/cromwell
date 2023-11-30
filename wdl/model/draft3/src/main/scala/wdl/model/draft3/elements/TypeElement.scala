package wdl.model.draft3.elements

import wom.types.WomPrimitiveType

sealed trait TypeElement extends LanguageElement

case class PrimitiveTypeElement(primitiveType: WomPrimitiveType) extends TypeElement

case class ArrayTypeElement(inner: TypeElement) extends TypeElement
case class MapTypeElement(keyType: TypeElement, valueType: TypeElement) extends TypeElement
case class OptionalTypeElement(maybeType: TypeElement) extends TypeElement
case class NonEmptyTypeElement(arrayType: TypeElement) extends TypeElement
case class PairTypeElement(leftType: TypeElement, rightType: TypeElement) extends TypeElement
case object ObjectTypeElement extends TypeElement

/**
  * Element to represent something like a Struct name which will eventually have to be linked back to a fixed WOM type.
  */
case class TypeAliasElement(alias: String) extends TypeElement
