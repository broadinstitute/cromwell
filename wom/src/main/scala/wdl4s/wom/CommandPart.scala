package wdl4s.wom

import wdl4s.wdl.values.WdlValue
import wdl4s.wom.expression.IoFunctionSet

trait CommandPart {
  def instantiate(inputsMap: Map[String, WdlValue],
                  functions: IoFunctionSet,
                  valueMapper: WdlValue => WdlValue): String
}
