package cromwell.binding

import cromwell.binding.types.WdlType

case class TaskOutput(name: String, wdlType: WdlType, expression: WdlExpression) {
  def evaluate: WdlValue = {
    val evaluation = expression.evaluate(???, ???)
    wdlType.checkCompatible(evaluation)
    evaluation
  }
}
