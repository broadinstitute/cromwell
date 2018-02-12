package wdl.model.draft3.elements

import common.validation.ErrorOr.ErrorOr
import simulacrum._

import scala.language.implicitConversions

case class WorkflowDefinitionElement(name: String,
                                     workflowOutputs: Seq[WorkflowOutputsElement]
                                    ) extends LanguageElement {
  override def children: Seq[LanguageElement] = workflowOutputs
}

@typeclass
trait WorkflowDefinitionElementMaker[A] {
  def convert(a: A): ErrorOr[WorkflowDefinitionElement]
}
