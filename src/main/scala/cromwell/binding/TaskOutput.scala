package cromwell.binding

import cromwell.binding.types.WdlType
import cromwell.binding.values.WdlValue

case class TaskOutput(name: String, wdlType: WdlType, expression: WdlExpression)
