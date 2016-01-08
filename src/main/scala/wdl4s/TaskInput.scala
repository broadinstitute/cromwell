package wdl4s

import wdl4s.types.WdlType

case class TaskInput(name: String, wdlType: WdlType, postfixQuantifier: Option[String] = None)
