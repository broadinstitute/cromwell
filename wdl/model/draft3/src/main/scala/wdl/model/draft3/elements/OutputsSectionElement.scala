package wdl.model.draft3.elements

final case class OutputsSectionElement(outputs: Seq[OutputDeclarationElement])
    extends WorkflowBodyElement
    with TaskSectionElement
