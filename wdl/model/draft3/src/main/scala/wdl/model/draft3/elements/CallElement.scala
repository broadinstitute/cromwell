package wdl.model.draft3.elements

final case class CallElement(callableName: String, alias: Option[String], body: Option[CallBodyElement])
  extends LanguageElement with WorkflowGraphElement
