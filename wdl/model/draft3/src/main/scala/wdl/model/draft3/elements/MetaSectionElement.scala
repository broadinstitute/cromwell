package wdl.model.draft3.elements

final case class MetaSectionElement(meta: Vector[Meta.KvPair]) extends WorkflowBodyElement with TaskSectionElement
