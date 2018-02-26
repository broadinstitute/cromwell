package wdl.model.draft3.elements

import wdl.model.draft3.elements.ExpressionElement.KvPair

final case class CallElement(callableName: String, alias: Option[String], body: Vector[KvPair]) extends LanguageElement with WorkflowGraphElement
