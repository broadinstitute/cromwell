package wdl4s.types
import wdl4s.values.WdlValue

/**
  * For when we need to assign a type to a composite wdl value that was created empty, may I present to you the least
  * (and yet at the same time, most) interesting of all the types, the WdlNothingType!
  */
case object WdlNothingType extends WdlType {
  override protected def coercion: PartialFunction[Any, WdlValue] = PartialFunction.empty
  override def toWdlString: String = "Nothing"
}
