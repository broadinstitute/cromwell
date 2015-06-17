package cromwell.binding

import cromwell.binding.types.WdlType

case class Declaration(wdlType: WdlType, name: String, expression: Option[WdlExpression])
