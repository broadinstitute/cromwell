package wdl.model.draft3.elements

import wom.SourceFileLocation

final case class ScatterElement(scatterName: String,
                                scatterExpression: ExpressionElement,
                                scatterVariableName: String,
                                graphElements: Seq[WorkflowGraphElement],
                                override val sourceLocation : Option[SourceFileLocation]) extends WorkflowGraphElement {

  // Scatter names do not contain intrinsic information about the scatter; rather they are a sort
  // of hash based on the declarations's physical location in the source.
  override def equals(other: scala.Any): Boolean = {
    other match {
      case otherScatter: ScatterElement =>
        this.scatterExpression == otherScatter.scatterExpression &&
        this.scatterVariableName == otherScatter.scatterVariableName &&
        this.graphElements == otherScatter.graphElements
      case _ => false
    }
  }

  // Shorthand to only include certain members for hashing purposes
  // https://stackoverflow.com/a/31915429/818054
  override def hashCode(): Int = (scatterExpression, scatterVariableName, graphElements).##
}
