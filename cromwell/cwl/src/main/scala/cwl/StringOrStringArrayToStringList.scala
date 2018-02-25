package cwl

import shapeless.Poly1

/**
  * From a logic standpoint it's easier to treat sources as a list of strings.
  */
object StringOrStringArrayToStringList extends Poly1{
  import Case._
  implicit def string: Aux[String, List[String]] = at[String]{List(_)}
  implicit def array: Aux[Array[String], List[String]] = at[Array[String]]{_.toList}
}
