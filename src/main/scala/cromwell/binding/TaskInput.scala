package cromwell.binding

import cromwell.binding.types.WdlType

case class TaskInput(name: String, types: Seq[WdlType], postfixQuantifier: Option[String])
