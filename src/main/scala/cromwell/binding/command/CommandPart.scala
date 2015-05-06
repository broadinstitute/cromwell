package cromwell.binding.command

import cromwell.binding.WdlValue

trait CommandPart {
  def instantiate(parameters: Map[String, WdlValue]): String
}
