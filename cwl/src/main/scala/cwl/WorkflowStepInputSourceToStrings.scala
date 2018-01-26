package cwl

import shapeless.Poly1

object WorkflowStepInputSourceToStrings extends Poly1{
  import Case._
  implicit def string: Aux[String, List[String]] = at[String]{List(_)}
  implicit def array: Aux[Array[String], List[String]] = at[Array[String]]{_.toList}
}
