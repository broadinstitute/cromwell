package wdl.model.draft3.elements

case class WorkflowDefinitionElement(name: String,
                                     inputsSection: Option[InputsSectionElement],
                                     outputsSection: Vector[OutputsSectionElement]
                                    ) extends LanguageElement
