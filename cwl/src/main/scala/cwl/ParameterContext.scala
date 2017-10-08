package cwl

import wdl.types._
import wdl.values._

object ParameterContext {
  val Empty = ParameterContext(
    inputs = WdlOptionalValue(WdlNothingType, None),
    self = WdlOptionalValue(WdlNothingType, None),
    runtime = WdlOptionalValue(WdlNothingType, None)
  )
}

case class ParameterContext(inputs: WdlValue, self: WdlValue, runtime: WdlValue)
