package wdl.model.draft3.elements

case class InputDeclarationElement(typeElement: TypeElement, name: String, expression: Option[Any]) extends LanguageElement
