package wdl.model.draft3.elements

final case class CallElement(callableReference: String, alias: Option[String], after: Option[String], body: Option[CallBodyElement])
  extends LanguageElement with WorkflowGraphElement
