package cromwell.binding

import cromwell.binding.command.ParameterCommandPart
import cromwell.binding.types.WdlType

case class WorkflowInput(fqn: FullyQualifiedName, wdlType: WdlType, postfixQuantifier: Option[String]) {
  val optional = postfixQuantifier match {
    case Some(s) => ParameterCommandPart.OptionalPostfixQuantifiers.contains(s)
    case _ => false
  }
}
