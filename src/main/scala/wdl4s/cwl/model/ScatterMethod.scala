package broad.cwl.model

import enumeratum._

sealed abstract class ScatterMethod(override val entryName: String) extends EnumEntry

object ScatterMethod extends Enum[ScatterMethod] {
  val values = findValues

  case object DotProduct extends ScatterMethod("dotproduct")
  case object NestedCrossProduct extends ScatterMethod("nested_crossproduct")
  case object FlatCrossProduct extends ScatterMethod("flat_crossproduct")
}
