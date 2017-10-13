package cwl

import wom.expression.IoFunctionSet
import wom.types.{WdlMapType, WdlNothingType, WdlStringType}
import wom.values.{WdlMap, WdlOptionalValue, WdlString, WdlValue}

object ParameterContext {
  val Empty = ParameterContext(
    inputs = WdlOptionalValue(WdlNothingType, None),
    self = WdlOptionalValue(WdlNothingType, None),
    runtime = WdlOptionalValue(WdlNothingType, None)
  )
}

case class ParameterContext(inputs: WdlValue, self: WdlValue, runtime: WdlValue) {
  def withInputs(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ParameterContext = {
    val wdlValueType = inputValues.values.headOption.map(_.wdlType).getOrElse(WdlNothingType)
    copy(
      inputs = WdlMap(
        WdlMapType(WdlStringType, wdlValueType),
        // TODO: WOM: convert inputValues (including WdlFile?) to inputs using the ioFunctionSet
        inputValues map { case (name, wdlValue) => WdlString(name) -> wdlValue }
      )
    )
  }
}
