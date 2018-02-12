package wdl.model.draft3.elements

import common.validation.ErrorOr.ErrorOr
import simulacrum.typeclass

case class WorkflowOutputDeclarationElement(name: String,
                                            valueType: String,
                                            value: String
                                ) extends LanguageElement {
  override def children: Seq[String] = Seq(name, valueType, value)
}

@typeclass
trait WorkflowOutputDeclarationElementMaker[A] {
  def convert(a: A): ErrorOr[WorkflowOutputDeclarationElement]
}
