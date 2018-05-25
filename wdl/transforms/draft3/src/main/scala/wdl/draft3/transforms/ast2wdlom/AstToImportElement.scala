package wdl.draft3.transforms.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import cats.instances.either._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.{Ast, AstNode}
import wdl.model.draft3.elements.{ImportElement, LanguageElement}
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._

object AstToImportElement {

  def convert(a: Ast): ErrorOr[ImportElement] = {
    val importPath: ErrorOr[String] = a.getAttributeAs[String]("uri").toValidated
    val alias: ErrorOr[Option[String]] = a.getAttributeAsOptional[String]("namespace").toValidated

    val aliasElementMaker: CheckedAtoB[AstNode, ImportStructRenameElement] = astNodeToAst andThen CheckedAtoB.fromErrorOr(convertAliasElement)
    val structRenames: ErrorOr[Vector[ImportStructRenameElement]] = a.getAttributeAsVector[ImportStructRenameElement]("aliases")(aliasElementMaker).toValidated
    val structRenameMap: ErrorOr[Map[String, String]] = structRenames.map(_.map(rename => rename.oldName -> rename.newName).toMap)

    (importPath, alias, structRenameMap) mapN { ImportElement }
  }

  private def convertAliasElement(a: Ast): ErrorOr[ImportStructRenameElement] = {
    val oldName: ErrorOr[String] = a.getAttributeAs[String]("old_name").toValidated
    val newName: ErrorOr[String] = a.getAttributeAs[String]("new_name").toValidated

    (oldName, newName) mapN ImportStructRenameElement
  }

  private final case class ImportStructRenameElement(oldName: String, newName: String) extends LanguageElement
}
