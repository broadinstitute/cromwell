package wdl.model.draft3.elements

final case class OutputDeclarationElement(typeElement: TypeElement, name: String, expression: String) extends WorkflowBodyElement with TaskBodyElement
