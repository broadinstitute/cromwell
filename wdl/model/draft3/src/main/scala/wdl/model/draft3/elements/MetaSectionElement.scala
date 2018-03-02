package wdl.model.draft3.elements

final case class MetaSectionElement(meta: Vector[MetaKvPair]) extends WorkflowBodyElement with TaskSectionElement
