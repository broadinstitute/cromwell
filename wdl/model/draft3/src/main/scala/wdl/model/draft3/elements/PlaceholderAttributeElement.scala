package wdl.model.draft3.elements

sealed trait PlaceholderAttributeElement extends LanguageElement

final case class DefaultAttributeElement(value: String) extends PlaceholderAttributeElement
final case class TrueAttributeElement(value: String) extends PlaceholderAttributeElement
final case class FalseAttributeElement(value: String) extends PlaceholderAttributeElement
final case class SepAttributeElement(value: String) extends PlaceholderAttributeElement

final case class PlaceholderAttributeSet(defaultAttribute: Option[String],
                                         trueAttribute: Option[String],
                                         falseAttribute: Option[String],
                                         sepAttribute: Option[String]
)

object PlaceholderAttributeSet {
  val empty = PlaceholderAttributeSet(None, None, None, None)
}
