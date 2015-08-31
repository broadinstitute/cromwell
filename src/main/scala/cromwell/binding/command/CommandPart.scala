package cromwell.binding.command

import cromwell.binding.expression.WdlFunctions
import cromwell.binding.values.WdlValue
import cromwell.binding.Declaration

trait CommandPart {
  def instantiate(declarations: Seq[Declaration], parameters: Map[String, WdlValue], functions: WdlFunctions[WdlValue]): String
}
