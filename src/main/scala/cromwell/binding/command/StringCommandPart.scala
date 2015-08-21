package cromwell.binding.command

import cromwell.binding.values.WdlValue
import cromwell.binding.{Declaration, WdlFunctions}

case class StringCommandPart(literal: String) extends CommandPart {
  override def toString: String = literal

  override def instantiate(declarations: Seq[Declaration], parameters: Map[String, WdlValue], functions: WdlFunctions): String = literal
}
