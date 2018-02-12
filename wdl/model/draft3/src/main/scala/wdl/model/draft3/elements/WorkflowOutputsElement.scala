package wdl.model.draft3.elements

import common.validation.ErrorOr.ErrorOr
import simulacrum.typeclass

class WorkflowOutputsElement(outputs: Seq[WorkflowOutputsElement]
                                       ) extends LanguageElement {
  override def children: Seq[LanguageElement] = outputs
}

@typeclass
trait WorkflowOutputsElementMaker[A] {
  def convert(a: A): ErrorOr[WorkflowOutputsElement]
}