package wdl.model.draft3.elements

final case class IfElement(conditionExpression: ExpressionElement,
                           graphElements: Seq[WorkflowGraphElement]) extends WorkflowGraphElement
