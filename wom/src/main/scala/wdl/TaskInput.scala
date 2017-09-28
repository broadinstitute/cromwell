package wdl

import wdl.types.WdlType

case class TaskInput(name: String, wdlType: WdlType, postfixQuantifier: Option[String] = None)
