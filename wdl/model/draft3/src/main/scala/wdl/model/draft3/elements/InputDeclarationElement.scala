package wdl.model.draft3.elements

case class InputDeclarationElement(typeElement: TypeElement, name: String, expression: Option[ExpressionElement]) extends LanguageElement
