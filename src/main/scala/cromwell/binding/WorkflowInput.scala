package cromwell.binding

import cromwell.binding.types.WdlType

case class WorkflowInput(fqn: FullyQualifiedName, types: Seq[WdlType], postfixQuantifier: Option[String]) {
  val optional = postfixQuantifier match {
    case Some(s) => Set("?", "*").contains(s)
    case _ => false
  }
}
