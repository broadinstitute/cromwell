package cwl

import wom.expression.IoFunctionSet
import wom.types.{WomMapType, WomNothingType, WomStringType}
import wom.values.{WomMap, WomOptionalValue, WomSingleFile, WomString, WomValue}

object ParameterContext {
  val Empty = ParameterContext()
}

case class ParameterContext(inputs: WomValue = WomOptionalValue(WomNothingType, None),
                            self: WomValue = WomOptionalValue(WomNothingType, None),
                            runtime: WomValue = WomOptionalValue(WomNothingType, None)) {
  def withInputs(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ParameterContext = {
    val womValueType = WomStringType
    copy(
      inputs = WomMap(
        WomMapType(WomStringType, womValueType),
        // TODO: WOM: convert inputValues (including WomFile?) to inputs using the ioFunctionSet
        inputValues map {
          case (name, WomSingleFile(path)) => WomString(name) -> WomString(path)
          case (name, womValue) => WomString(name) -> womValue
        }
      )
    )
  }
}
