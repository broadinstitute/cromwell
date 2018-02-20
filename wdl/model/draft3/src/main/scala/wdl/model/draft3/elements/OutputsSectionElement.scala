package wdl.model.draft3.elements

final case class OutputsSectionElement(outputs: Seq[SetDeclarationElement]) extends WorkflowBodyElement with TaskBodyElement
