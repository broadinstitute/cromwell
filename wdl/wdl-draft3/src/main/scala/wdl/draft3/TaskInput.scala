package wdl.draft3

import wom.types.WomType

case class TaskInput(name: String, womType: WomType, postfixQuantifier: Option[String] = None)
