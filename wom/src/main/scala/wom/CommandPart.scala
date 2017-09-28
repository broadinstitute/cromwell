package wom

import wdl.values.WdlValue
import wom.expression.IoFunctionSet

trait CommandPart {
  def instantiate(inputsMap: Map[String, WdlValue],
                  functions: IoFunctionSet,
                  valueMapper: WdlValue => WdlValue): String
}
