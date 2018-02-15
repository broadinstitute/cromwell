package wdl.model.draft3.elements

final case class WorkflowDefinitionElement(name: String,
                                           outputsSection: Vector[OutputsSectionElement]
                                    ) extends LanguageElement {
  override def children: Seq[LanguageElement] = outputsSection
}
