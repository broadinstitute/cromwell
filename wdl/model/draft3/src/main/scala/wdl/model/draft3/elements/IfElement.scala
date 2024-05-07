package wdl.model.draft3.elements

final case class IfElement(conditionExpression: ExpressionElement, graphElements: Seq[WorkflowGraphElement])
    extends WorkflowGraphElement {
  override def toString: String = s"""Condition "$conditionExpression""""
}
