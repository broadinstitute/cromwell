package wdl4s.command

import wdl4s.expression.WdlFunctions
import wdl4s.values.WdlValue
import wdl4s.Declaration

case class StringCommandPart(literal: String) extends CommandPart {
  override def toString: String = literal

  override def instantiate(declarations: Seq[Declaration],
                           parameters: Map[String, WdlValue],
                           functions: WdlFunctions[WdlValue],
                           valueMapper: WdlValue => WdlValue = (v) => v): String = literal
}
