package wdl.model.draft3.elements

final case class WorkflowOutputDeclarationElement(valueType: String, name: String, value: String) extends LanguageElement {
  override def children: Seq[LanguageElement] = Seq.empty
}
