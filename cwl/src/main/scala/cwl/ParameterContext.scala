package cwl

import wom.expression.IoFunctionSet
import wom.types.{WomMapType, WomNothingType, WomStringType}
import wom.values.{WomMap, WomOptionalValue, WomString, WomValue}

object ParameterContext {
  val Empty = ParameterContext(
    inputs = WomOptionalValue(WomNothingType, None),
    self = WomOptionalValue(WomNothingType, None),
    runtime = WomOptionalValue(WomNothingType, None)
  )
}

case class ParameterContext(inputs: WomValue, self: WomValue, runtime: WomValue) {
  def withInputs(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ParameterContext = {
    val wdlValueType = inputValues.values.headOption.map(_.womType).getOrElse(WomNothingType)
    copy(
      inputs = WomMap(
        WomMapType(WomStringType, wdlValueType),
        // TODO: WOM: convert inputValues (including WdlFile?) to inputs using the ioFunctionSet
        inputValues map { case (name, womValue) => WomString(name) -> womValue }
      )
    )
  }
}
