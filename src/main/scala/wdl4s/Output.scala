package wdl4s

trait Output extends DeclarationInterface {
  def requiredExpression: WdlExpression
  
  override val postfixQuantifier = None
  override val expression = Option(requiredExpression)
}
