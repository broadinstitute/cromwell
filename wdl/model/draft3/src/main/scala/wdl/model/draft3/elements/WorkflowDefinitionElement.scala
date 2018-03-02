package wdl.model.draft3.elements

final case class WorkflowDefinitionElement(name: String,
                                           inputsSection: Option[InputsSectionElement],
                                           graphElements: Set[WorkflowGraphElement],
                                           outputsSection: Option[OutputsSectionElement],
                                           meta: Option[MetaSectionElement],
                                           parameterMeta: Option[ParameterMetaSectionElement]) extends FileBodyElement
