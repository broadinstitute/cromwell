package wdl.model.draft3.elements

final case class ParameterMetaSectionElement(metaAttributes: Map[String, MetaValueElement]) extends WorkflowBodyElement with TaskSectionElement
