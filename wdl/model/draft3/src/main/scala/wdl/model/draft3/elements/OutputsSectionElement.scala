package wdl.model.draft3.elements

final case class OutputsSectionElement(outputs: Seq[DeclarationElement]) extends WorkflowBodyElement with TaskBodyElement
