package wdl.command

import wdl._
import wdl.expression.WdlFunctions
import wom.values.WdlValue

case class StringCommandPart(literal: String) extends WdlCommandPart {
  override def toString: String = literal

  override def instantiate(declarations: Seq[Declaration], inputsMap: Map[String, WdlValue], functions: WdlFunctions[WdlValue], valueMapper: (WdlValue) => WdlValue): String = literal
}
