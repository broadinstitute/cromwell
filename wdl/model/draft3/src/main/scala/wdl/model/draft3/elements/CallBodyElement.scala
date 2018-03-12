package wdl.model.draft3.elements

import wdl.model.draft3.elements.ExpressionElement.KvPair

final case class CallBodyElement(inputs: Vector[KvPair]) extends LanguageElement
