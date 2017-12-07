package cwl

import shapeless.Poly1

object ScatterLogic {

  object ScatterVariablesPoly extends Poly1 {
    implicit def fromString: Case.Aux[String, List[String]] = at[String] { s: String => List(s) }
    implicit def fromStringList: Case.Aux[Array[String], List[String]] = at[Array[String]] { l: Array[String] => l.toList }
  }
}
