package wom.types

import wom.values.WomValue

/**
  * For when we need to assign a type to a composite wom value that was created empty, may I present to you the least
  * (and yet at the same time, most) interesting of all the types, the WomNothingType!
  */
case object WomNothingType extends WomType {
  override def coercion(): PartialFunction[Any, WomValue] = PartialFunction.empty
  override def toDisplayString: String = "Nothing"
}
