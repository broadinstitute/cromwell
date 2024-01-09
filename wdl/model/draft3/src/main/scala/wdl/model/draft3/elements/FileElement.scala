package wdl.model.draft3.elements

final case class FileElement(imports: Seq[ImportElement],
                             structs: Seq[StructElement],
                             workflows: Seq[WorkflowDefinitionElement],
                             tasks: Seq[TaskDefinitionElement]
) extends LanguageElement
