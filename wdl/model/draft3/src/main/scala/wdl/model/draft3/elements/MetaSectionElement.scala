package wdl.model.draft3.elements

final case class MetaSectionElement(meta: Map[String, MetaValueElement]) extends WorkflowBodyElement with TaskSectionElement
