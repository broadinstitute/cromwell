package wdl4s

import wdl4s.types.WdlType

case class WorkflowInput(fqn: FullyQualifiedName, wdlType: WdlType, postfixQuantifier: Option[String]) {
  val optional = postfixQuantifier match {
    case Some(s) if s == "?" => true
    case _ => false
  }
}
