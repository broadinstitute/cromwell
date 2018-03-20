package cwl

import wom.callable.RuntimeEnvironment
import wom.types.WomNothingType
import wom.values.{WomOptionalValue, WomValue}

object ParameterContext {
  val EmptySelf = WomOptionalValue(WomNothingType, None)
}

case class ParameterContext(inputs: Map[String, WomValue] = Map.empty,
                            self: WomValue = ParameterContext.EmptySelf,
                            runtimeOption: Option[RuntimeEnvironment] = None)
