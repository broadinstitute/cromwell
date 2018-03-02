package wdl.model.draft3.elements

final case class TaskDefinitionElement(name: String,
                                       inputsSection: Option[InputsSectionElement],
                                       outputsSection: Option[OutputsSectionElement],
                                       commandSection: CommandSectionElement,
                                       runtimeSection: Option[RuntimeAttributesSectionElement],
                                       meta: Option[MetaSectionElement],
                                       parameterMeta: Option[ParameterMetaSectionElement]) extends FileBodyElement
