package wdl.model.draft3.elements

final case class WorkflowDefinitionElement(name: String,
                                     workflowOutputs: Vector[WorkflowOutputsElement]
                                    ) extends LanguageElement {
  override def children: Seq[LanguageElement] = workflowOutputs
}
