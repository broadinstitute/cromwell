package wdl.model.draft3.elements

final case class ScatterElement(scatterExpression: ExpressionElement,
                                scatterVariableName: String,
                                graphElements: Seq[WorkflowGraphElement]) extends WorkflowGraphElement
