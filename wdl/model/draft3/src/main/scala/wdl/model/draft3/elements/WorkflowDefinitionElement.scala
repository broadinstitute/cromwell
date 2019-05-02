package wdl.model.draft3.elements
import wom.LexicalInformation

final case class WorkflowDefinitionElement(name: String,
                                           inputsSection: Option[InputsSectionElement],
                                           graphElements: Set[WorkflowGraphElement],
                                           outputsSection: Option[OutputsSectionElement],
                                           metaSection: Option[MetaSectionElement],
                                           parameterMetaSection: Option[ParameterMetaSectionElement],
                                           override val lexInfo : Option[LexicalInformation] = None) extends FileBodyElement
