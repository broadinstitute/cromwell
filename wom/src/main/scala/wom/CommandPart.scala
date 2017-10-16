package wom

import wom.expression.IoFunctionSet
import wom.graph.LocalName
import wom.values.WomValue

trait CommandPart {
  def instantiate(inputsMap: Map[LocalName, WomValue],
                  functions: IoFunctionSet,
                  valueMapper: WomValue => WomValue): String
}
