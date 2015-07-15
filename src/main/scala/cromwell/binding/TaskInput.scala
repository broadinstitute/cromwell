package cromwell.binding

import cromwell.binding.types.WdlType

case class TaskInput(name: String, wdlType: WdlType, postfixQuantifier: Option[String] = None)
