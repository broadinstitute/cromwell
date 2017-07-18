package wdl4s.wdl.command

import wdl4s.wdl.expression.WdlFunctions
import wdl4s.wdl.values.WdlValue
import wdl4s.wdl._

trait CommandPart {
  def instantiate(declarations: Seq[Declaration],
                  inputsMap: EvaluatedTaskInputs,
                  functions: WdlFunctions[WdlValue],
                  valueMapper: WdlValue => WdlValue): String
}
