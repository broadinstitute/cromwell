package wdl4s.wdl.command

import wdl4s.wdl._
import wdl4s.wdl.expression.WdlFunctions
import wdl4s.wdl.values.WdlValue

case class StringCommandPart(literal: String) extends WdlCommandPart {
  override def toString: String = literal

  override def instantiate(declarations: Seq[Declaration], inputsMap: Map[String, WdlValue], functions: WdlFunctions[WdlValue], valueMapper: (WdlValue) => WdlValue): String = literal
}
