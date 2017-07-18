package wdl4s.cwl

object ScatterMethod extends Enumeration {
  type ScatterMethod = Value

  val DotProduct = Value("dotproduct")
  val NestedCrossProduct = Value("nested_crossproduct")
  val FlatCrossProduct = Value("flat_crossproduct")
}
