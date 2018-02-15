package wdl.model.draft3.elements

case class FileElement(imports: Seq[ImportElement],
                       workflows: Seq[WorkflowDefinitionElement],
                       tasks: Seq[TaskDefinitionElement]) extends LanguageElement
