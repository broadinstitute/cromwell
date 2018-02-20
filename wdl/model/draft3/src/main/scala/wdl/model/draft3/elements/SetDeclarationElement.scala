package wdl.model.draft3.elements

/**
  * A Declaration which has a setter (ie not a declaration of a set!)
  */
final case class SetDeclarationElement(typeElement: TypeElement, name: String, expression: ExpressionElement) extends TaskBodyElement with WorkflowGraphNodeElement
