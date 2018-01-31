package wdl.draft2

trait Output extends DeclarationInterface {
  def requiredExpression: WdlExpression

  override val expression = Option(requiredExpression)
}
