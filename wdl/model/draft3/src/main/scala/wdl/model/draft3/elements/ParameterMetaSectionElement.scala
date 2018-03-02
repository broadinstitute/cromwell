package wdl.model.draft3.elements

final case class ParameterMetaSectionElement(metaAttributes: Vector[Meta.KvPair]) extends WorkflowBodyElement with TaskSectionElement
