package wdl.draft2.model

import wom.types.WomType

case class TaskInput(name: String, womType: WomType, postfixQuantifier: Option[String] = None)
