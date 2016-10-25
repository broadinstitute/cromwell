package wdl4s.command

import wdl4s.expression.WdlFunctions
import wdl4s.values.WdlValue
import wdl4s._

trait CommandPart {
  def instantiate(declarations: Seq[Declaration],
                  inputsMap: EvaluatedTaskInputs,
                  functions: WdlFunctions[WdlValue],
                  valueMapper: WdlValue => WdlValue): String
}
