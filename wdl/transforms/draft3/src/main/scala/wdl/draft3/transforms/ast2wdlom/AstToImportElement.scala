package wdl.draft3.transforms.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._

import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.Ast
import wdl.model.draft3.elements.ImportElement
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._

object AstToImportElement {

  def convert(a: Ast): ErrorOr[ImportElement] = {
    val importPath: ErrorOr[String] = a.getAttributeAs[String]("uri").toValidated
    val alias: ErrorOr[Option[String]] = a.getAttributeAsOptional[String]("namespace").toValidated

    (importPath, alias) mapN { ImportElement }
  }
}
