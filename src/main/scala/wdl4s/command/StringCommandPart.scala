package wdl4s.command

import wdl4s._
import wdl4s.expression.WdlFunctions
import wdl4s.values.WdlValue

case class StringCommandPart(literal: String) extends CommandPart {
  override def toString: String = literal

  override def instantiate(declarations: Seq[Declaration], inputsMap: EvaluatedTaskInputs, functions: WdlFunctions[WdlValue], valueMapper: (WdlValue) => WdlValue): String = literal
}
