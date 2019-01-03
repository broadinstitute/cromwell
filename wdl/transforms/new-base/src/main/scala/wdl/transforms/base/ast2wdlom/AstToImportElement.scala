package wdl.transforms.base.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import cats.instances.either._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.{ImportElement, LanguageElement, StaticString}

object AstToImportElement {

  def astToImportElement(implicit astNodeToStaticString: CheckedAtoB[GenericAstNode, StaticString]): CheckedAtoB[GenericAst, ImportElement] =
    CheckedAtoB.fromErrorOr("convert Ast to ImportElement") { a =>

      val importPath: ErrorOr[String] = a.getAttributeAs[StaticString]("uri").map(_.value).toValidated
      val alias: ErrorOr[Option[String]] = a.getAttributeAsOptional[String]("namespace").toValidated

      val aliasElementMaker: CheckedAtoB[GenericAstNode, ImportStructRenameElement] = astNodeToAst andThen CheckedAtoB.fromErrorOr(convertAliasElement _)
      val structRenames: ErrorOr[Vector[ImportStructRenameElement]] = a.getAttributeAsVector[ImportStructRenameElement]("aliases")(aliasElementMaker).toValidated
      val structRenameMap: ErrorOr[Map[String, String]] = structRenames.map(_.map(rename => rename.oldName -> rename.newName).toMap)

      (importPath, alias, structRenameMap) mapN ImportElement
    }

  private def convertAliasElement(a: GenericAst): ErrorOr[ImportStructRenameElement] = {
    val oldName: ErrorOr[String] = a.getAttributeAs[String]("old_name").toValidated
    val newName: ErrorOr[String] = a.getAttributeAs[String]("new_name").toValidated

    (oldName, newName) mapN ImportStructRenameElement
  }

  private final case class ImportStructRenameElement(oldName: String, newName: String) extends LanguageElement
}
