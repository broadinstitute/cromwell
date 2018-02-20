package wdl.model.draft3.elements

/**
  * A Declaration which has a setter
  */
final case class DeclarationElement(typeElement: TypeElement, name: String, expression: ExpressionElement) extends TaskBodyElement with WorkflowGraphNodeElement
