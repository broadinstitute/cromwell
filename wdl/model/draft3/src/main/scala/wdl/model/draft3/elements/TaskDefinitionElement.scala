package wdl.model.draft3.elements

final case class TaskDefinitionElement(name: String, inputsSection: Option[InputsSectionElement],
                                       outputsSection: Option[OutputsSectionElement],
                                       commandSection: CommandSectionElement,
                                       runtimeSection: Option[RuntimeAttributesSectionElement],
                                       meta: Map[String, MetaValueElement],
                                       parameterMeta: Map[String, MetaValueElement]) extends FileBodyElement
