package wdl.model.draft3.elements

case class WorkflowDefinitionElement(name: String,
                                     inputsSection: Option[InputsSectionElement],
                                     graphElements: Set[WorkflowGraphElement],
                                     outputsSection: Option[OutputsSectionElement]
                                    ) extends FileBodyElement
