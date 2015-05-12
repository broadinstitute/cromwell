package cromwell.binding.command

import cromwell.binding.values.WdlValue

case class StringCommandPart(literal: String) extends CommandPart {
  override def toString: String = literal

  def instantiate(parameters: Map[String, WdlValue]): String = literal
}
