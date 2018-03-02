package wdl.model.draft3.elements

final case class ParameterMetaSectionElement(metaAttributes: Vector[MetaKvPair]) extends WorkflowBodyElement with TaskSectionElement
