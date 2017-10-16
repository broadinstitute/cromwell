package wdl

import wom.types.WomType

case class TaskInput(name: String, womType: WomType, postfixQuantifier: Option[String] = None)
