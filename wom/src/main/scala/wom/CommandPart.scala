package wom

import wom.expression.IoFunctionSet
import wom.graph.LocalName
import wom.values.WdlValue

trait CommandPart {
  def instantiate(inputsMap: Map[LocalName, WdlValue],
                  functions: IoFunctionSet,
                  valueMapper: WdlValue => WdlValue): String
}
