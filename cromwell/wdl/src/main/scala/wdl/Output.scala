package wdl

trait Output extends DeclarationInterface {
  def requiredExpression: WdlExpression

  override val expression = Option(requiredExpression)
}
