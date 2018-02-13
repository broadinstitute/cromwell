package wdl.model.draft3.elements

import common.validation.ErrorOr.ErrorOr
import simulacrum._

import scala.language.implicitConversions

case class ImportElement(importUrl: String,
                         alias: Option[String]) extends LanguageElement {
  override def children: Seq[LanguageElement] = List.empty
}

@typeclass
trait ImportElementMaker[A] {
  def convert(a: A): ErrorOr[ImportElement]
}
