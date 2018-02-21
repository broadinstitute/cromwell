package wdl.model.draft3.elements

case class WorkflowDefinitionElement(name: String,
                                     inputsSection: Option[InputsSectionElement],
                                     graphElements: Set[WorkflowGraphNodeElement],
                                     outputsSection: Option[OutputsSectionElement]
                                    ) extends FileBodyElement
