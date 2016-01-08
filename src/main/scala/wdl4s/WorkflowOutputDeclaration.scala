package wdl4s

case class WorkflowOutputDeclaration(fqn: String, wildcard: Boolean) {

  def outputMatchesDeclaration(outputFqn: String, wildcardsAllowed: Boolean): Boolean = {
    if (wildcard) {
      wildcardsAllowed && outputFqn.startsWith(fqn)
    } else {
      outputFqn.equals(fqn)
    }
  }
}
