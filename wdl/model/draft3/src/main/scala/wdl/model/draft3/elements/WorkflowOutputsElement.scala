package wdl.model.draft3.elements

final case class WorkflowOutputsElement(outputs: Seq[WorkflowOutputDeclarationElement]) extends LanguageElement {
  override def children: Seq[LanguageElement] = outputs
}
