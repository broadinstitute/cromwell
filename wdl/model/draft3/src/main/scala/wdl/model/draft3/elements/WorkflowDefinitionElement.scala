package wdl.model.draft3.elements

import common.validation.ErrorOr.ErrorOr
import simulacrum._

import scala.language.implicitConversions

case class WorkflowDefinitionElement(identifier: String
                                    // TODO: sections
                                    ) extends LanguageElement {
  override def children: Seq[LanguageElement] = List.empty // TODO: sections
}

@typeclass
trait WorkflowDefinitionElementMaker[A] {
  def convert(a: A): ErrorOr[WorkflowDefinitionElement]
}
