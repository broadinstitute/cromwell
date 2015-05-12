package cromwell.binding.command

import cromwell.binding.values.WdlValue

trait CommandPart {
  def instantiate(parameters: Map[String, WdlValue]): String
}
