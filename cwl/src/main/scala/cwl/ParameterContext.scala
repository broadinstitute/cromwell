package cwl

import wom.expression.IoFunctionSet
import wom.types.{WomMapType, WomNothingType, WomStringType}
import wom.values.{WomMap, WomOptionalValue, WomSingleFile, WomString, WomValue}

object ParameterContext {
  val Empty = ParameterContext(
    inputs = WomOptionalValue(WomNothingType, None),
    self = WomOptionalValue(WomNothingType, None),
    runtime = WomOptionalValue(WomNothingType, None)
  )
}

case class ParameterContext(inputs: WomValue, self: WomValue, runtime: WomValue) {
  def withInputs(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ParameterContext = {
    val wdlValueType = WomStringType
    copy(
      inputs = WomMap(
        WomMapType(WomStringType, wdlValueType),
        // TODO: WOM: convert inputValues (including WdlFile?) to inputs using the ioFunctionSet
        inputValues map {
          case (name, WomSingleFile(path)) => WomString(name) -> WomString(path)
          case (name, womValue) => WomString(name) -> womValue
        }
      )
    )
  }
}
