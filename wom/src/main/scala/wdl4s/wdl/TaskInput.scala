package wdl4s.wdl

import wdl4s.wdl.types.WdlType

case class TaskInput(name: String, wdlType: WdlType, postfixQuantifier: Option[String] = None)
