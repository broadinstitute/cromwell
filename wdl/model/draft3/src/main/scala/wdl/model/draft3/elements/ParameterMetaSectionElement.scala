package wdl.model.draft3.elements

import wom.callable.MetaValueElement

final case class ParameterMetaSectionElement(metaAttributes: Map[String, MetaValueElement]) extends WorkflowBodyElement with TaskSectionElement
