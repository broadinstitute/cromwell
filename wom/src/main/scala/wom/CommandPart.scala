package wom

import wom.callable.RuntimeEnvironment
import wom.expression.IoFunctionSet
import wom.graph.LocalName
import wom.values.WomValue

trait CommandPart {
  def instantiate(inputsMap: Map[LocalName, WomValue],
                  functions: IoFunctionSet,
                  valueMapper: WomValue => WomValue,
                  runtimeEnvironment: RuntimeEnvironment): String
}
