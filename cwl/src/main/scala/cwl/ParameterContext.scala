package cwl

import wom.callable.RuntimeEnvironment
import wom.types.WomNothingType
import wom.values.{WomOptionalValue, WomValue}

case class ParameterContext(inputs: Map[String, WomValue] = Map.empty,
                            self: WomValue = WomOptionalValue(WomNothingType, None),
                            runtimeOption: Option[RuntimeEnvironment] = None)
