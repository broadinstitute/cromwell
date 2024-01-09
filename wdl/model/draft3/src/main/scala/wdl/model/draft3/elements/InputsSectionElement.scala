package wdl.model.draft3.elements

final case class InputsSectionElement(inputDeclarations: Seq[InputDeclarationElement])
    extends WorkflowBodyElement
    with TaskSectionElement
