package wdl

import wom.types.WdlType

case class TaskInput(name: String, wdlType: WdlType, postfixQuantifier: Option[String] = None)
