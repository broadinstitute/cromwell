package wdl.transforms.base.ast2wdlom

import cats.syntax.either._
import cats.syntax.apply._
import cats.instances.either._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.{StructElement, StructEntryElement, TypeElement}

object AstToStructElement {
  def astToStructElement(implicit astNodeToTypeElement: CheckedAtoB[GenericAstNode, TypeElement]
                         ): CheckedAtoB[GenericAst, StructElement] = CheckedAtoB.fromErrorOr("convert AST to struct definition") { a =>

    def convertAstToStructEntry(a: GenericAst): ErrorOr[StructEntryElement] = {
      val name: ErrorOr[String] = a.getAttributeAs[String]("name").toValidated
      val typeElement: ErrorOr[TypeElement] = a.getAttributeAs[TypeElement]("type").toValidated

      (name, typeElement).mapN(StructEntryElement.apply)
    }

    implicit val astNodeToStructEntry: CheckedAtoB[GenericAstNode, StructEntryElement] = astNodeToAst andThen CheckedAtoB.fromErrorOr("convert AST to struct entry")(convertAstToStructEntry)
    val nameValidation: ErrorOr[String] = a.getAttributeAs[String]("name").toValidated
    val entriesValidation: ErrorOr[Vector[StructEntryElement]] = a.getAttributeAsVector[StructEntryElement]("entries").toValidated

    (nameValidation, entriesValidation) mapN { (name, entries) => StructElement(name, entries) }


  }
}
