package wdl.model.draft3.elements

final case class TaskDefinitionElement(name: String,
                                       inputsSection: Option[InputsSectionElement],
                                       outputsSection: Option[OutputsSectionElement],
                                       commandSection: CommandSectionElement,
                                       runtimeSection: Option[RuntimeAttributesSectionElement],
                                       metaSection: Option[MetaSectionElement],
                                       parameterMetaSection: Option[ParameterMetaSectionElement]) extends FileBodyElement
