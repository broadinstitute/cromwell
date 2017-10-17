package cwl

import wdl.types._
import wdl.values._
import wom.expression.IoFunctionSet

object ParameterContext {
  val Empty = ParameterContext(
    inputs = WdlOptionalValue(WdlNothingType, None),
    self = WdlOptionalValue(WdlNothingType, None),
    runtime = WdlOptionalValue(WdlNothingType, None)
  )
}

case class ParameterContext(inputs: WdlValue, self: WdlValue, runtime: WdlValue) {
  def withInputs(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ParameterContext = {
    val wdlValueType = WdlStringType
    copy(
      inputs = WdlMap(
        WdlMapType(WdlStringType, wdlValueType),
        // TODO: WOM: convert inputValues (including WdlFile?) to inputs using the ioFunctionSet
        inputValues map {
          case (name, WdlSingleFile(path)) => WdlString(name) -> WdlString(path)
          case (name, wdlValue) => WdlString(name) -> wdlValue
        }
      )
    )
  }
}
