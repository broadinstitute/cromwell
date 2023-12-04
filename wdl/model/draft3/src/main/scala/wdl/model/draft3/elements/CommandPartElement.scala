package wdl.model.draft3.elements

sealed trait CommandPartElement extends TaskSectionElement

object CommandPartElement {
  final case class StringCommandPartElement(value: String) extends CommandPartElement
  final case class PlaceholderCommandPartElement(expressionElement: ExpressionElement,
                                                 attributes: PlaceholderAttributeSet
  ) extends CommandPartElement
}
