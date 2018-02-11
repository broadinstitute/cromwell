package wdl.draft2.model

import wdl.draft2.parser.WdlParser.Ast

case class WorkflowOutputWildcard(fqn: String, wildcard: Boolean, ast: Ast) {

  def outputMatchesDeclaration(outputFqn: String, wildcardsAllowed: Boolean): Boolean = {
    if (wildcard) {
      val callFqn = outputFqn.substring(0, outputFqn.lastIndexOf('.'))
      wildcardsAllowed && callFqn.equals(fqn)
    } else {
      outputFqn.equals(fqn)
    }
  }
}
