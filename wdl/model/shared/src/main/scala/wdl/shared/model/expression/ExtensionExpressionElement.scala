package wdl.shared.model.expression

/**
 * An expression element that is not part of the shared WDL 1.0+ element hierarchy, e.g. an operator or function added
 * by a newer WDL version. It describes its own sub-expressions, so that version-independent code (such as identifier
 * renaming and writing WDL source) can process it without knowing its concrete type.
 *
 * @tparam E the expression element type
 */
trait ExtensionExpressionElement[E] {

  /** The sub-expressions of this element, in order. */
  def params: Seq[E]

  /** A copy of this element with its sub-expressions replaced, in the same order as `params`. */
  def withParams(newParams: Seq[E]): E

  /** This element as WDL source, given a way to write its sub-expressions. */
  def toWdlV1(paramToWdl: E => String): String
}
