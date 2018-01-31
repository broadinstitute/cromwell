package wdl.draft3

trait Output extends DeclarationInterface {
  def requiredExpression: WdlExpression

  override val expression = Option(requiredExpression)
}
