package cromwell.binding.command

import cromwell.binding.values.WdlValue
import cromwell.binding.{Declaration, WdlFunctions}

trait CommandPart {
  def instantiate(declarations: Seq[Declaration], parameters: Map[String, WdlValue], functions: WdlFunctions): String
}
