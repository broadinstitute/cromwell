package wdl.model.draft3.elements

case class WorkflowDefinitionElement(name: String,
                                     inputsSection: Option[InputsSectionElement],
                                     outputsSection: Option[OutputsSectionElement]
                                    ) extends LanguageElement
