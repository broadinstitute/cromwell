package wdl.command

import wdl._
import wdl.expression.WdlFunctions
import wom.values.WomValue

case class StringCommandPart(literal: String) extends WdlCommandPart {
  override def toString: String = literal

  override def instantiate(declarations: Seq[Declaration], inputsMap: Map[String, WomValue], functions: WdlFunctions[WomValue], valueMapper: (WomValue) => WomValue): String = literal
}
