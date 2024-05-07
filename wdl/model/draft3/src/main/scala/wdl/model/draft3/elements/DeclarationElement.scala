package wdl.model.draft3.elements

/**
  * Content of an intermediate or output declaration:
  */
final case class DeclarationContent(typeElement: TypeElement, name: String, expression: ExpressionElement)

/**
  * A Declaration outside of an input or output block
  */
final case class IntermediateValueDeclarationElement(typeElement: TypeElement,
                                                     name: String,
                                                     expression: ExpressionElement
) extends WorkflowGraphElement
    with TaskSectionElement {
  override def toString: String = s"""Declaration "$name" ($typeElement)"""
}

object IntermediateValueDeclarationElement {
  def fromContent(content: DeclarationContent): IntermediateValueDeclarationElement =
    IntermediateValueDeclarationElement(content.typeElement, content.name, content.expression)
}

/**
  * A declaration in an output block
  */
final case class OutputDeclarationElement(typeElement: TypeElement, name: String, expression: ExpressionElement)
    extends LanguageElement
    with WorkflowGraphElement {
  override def toString: String = s"""Output "$name" ($typeElement)"""
}

object OutputDeclarationElement {
  def fromContent(content: DeclarationContent): OutputDeclarationElement =
    OutputDeclarationElement(content.typeElement, content.name, content.expression)
}

/**
  * A declaration in an input block
  */
final case class InputDeclarationElement(typeElement: TypeElement, name: String, expression: Option[ExpressionElement])
    extends LanguageElement
    with WorkflowGraphElement {
  override def toString: String = s"""Input "$name" ($typeElement)"""
}

object DeclarationElement {
  /* Custom unapply so that elsewhere we can do things like this, and otherwise treat all declarations the same:
  languageElement match {
    case DeclarationElement(typeElement, name, Some(expr)) => ...
  }
   */
  def unapply(languageElement: LanguageElement): Option[(TypeElement, String, Option[ExpressionElement])] =
    languageElement match {
      case IntermediateValueDeclarationElement(typeElement, name, expr) => Option((typeElement, name, Option(expr)))
      case OutputDeclarationElement(typeElement, name, expr) => Option((typeElement, name, Option(expr)))
      case InputDeclarationElement(typeElement, name, expr) => Option((typeElement, name, expr))
      case _ => None
    }
}
