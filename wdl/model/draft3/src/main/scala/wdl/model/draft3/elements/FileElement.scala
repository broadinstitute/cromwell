package wdl.model.draft3.elements

import common.validation.ErrorOr.ErrorOr
import simulacrum._

import scala.language.implicitConversions

case class FileElement(imports: Seq[ImportElement],
                       workflows: Seq[WorkflowDefinitionElement],
                       tasks: Seq[TaskDefinitionElement]) extends LanguageElement {

  override def children: Seq[LanguageElement] = imports ++ workflows ++ tasks
}

@typeclass
trait FileElementMaker[A] {
  def convert(a: A): ErrorOr[FileElement]
}
