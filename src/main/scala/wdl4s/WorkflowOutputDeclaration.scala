package wdl4s

case class WorkflowOutputDeclaration(fqn: String, wildcard: Boolean) {

  def outputMatchesDeclaration(outputFqn: String, wildcardsAllowed: Boolean): Boolean = {
    if (wildcard) {
      val callFqn = outputFqn.substring(0, outputFqn.lastIndexOf('.'))
      wildcardsAllowed && callFqn.equals(fqn)
    } else {
      outputFqn.equals(fqn)
    }
  }
}
