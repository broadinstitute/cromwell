package wdl4s

import wdl4s.types.WdlType

case class WorkflowInput(fqn: FullyQualifiedName, wdlType: WdlType, postfixQuantifier: Option[String]) {
  val optional = postfixQuantifier match {
    case Some(s) if s == Declaration.OptionalPostfixQuantifier => true
    case _ => false
  }
}
