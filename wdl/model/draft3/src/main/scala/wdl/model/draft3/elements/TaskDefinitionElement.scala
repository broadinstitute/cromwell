package wdl.model.draft3.elements
import wom.SourceFileLocation

final case class TaskDefinitionElement(name: String,
                                       inputsSection: Option[InputsSectionElement],
                                       declarations: Seq[IntermediateValueDeclarationElement],
                                       outputsSection: Option[OutputsSectionElement],
                                       commandSection: CommandSectionElement,
                                       runtimeSection: Option[RuntimeAttributesSectionElement],
                                       metaSection: Option[MetaSectionElement],
                                       parameterMetaSection: Option[ParameterMetaSectionElement],
                                       override val sourceLocation : Option[SourceFileLocation]) extends FileBodyElement
