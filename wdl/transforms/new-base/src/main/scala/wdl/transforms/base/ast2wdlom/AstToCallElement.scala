package wdl.transforms.base.ast2wdlom

import cats.syntax.validated._
import cats.syntax.either._
import cats.instances.either._
import cats.syntax.apply._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.{CallBodyElement, CallElement}
import wdl.model.draft3.elements.ExpressionElement.KvPair

object AstToCallElement {

  def astToCallElement(implicit astNodeToKvPair: CheckedAtoB[GenericAstNode, KvPair]): CheckedAtoB[GenericAst, CallElement] = CheckedAtoB.fromErrorOr { ast =>
    def convertBodyElement(a: GenericAst): Checked[CallBodyElement] = {
      a.getAttributeAsVector[KvPair]("inputs") map CallBodyElement
    }

    val callableNameValidation: ErrorOr[String] = astNodeToString(ast.getAttribute("task")).toValidated

    val aliasValidation: ErrorOr[Option[String]] = Option(ast.getAttribute("alias")) match {
      case Some(a) => astNodeToString(a).map(Some.apply).toValidated
      case None => None.validNel
    }

    val afterValidation: ErrorOr[Option[String]] = Option(ast.getAttribute("after")) match {
      case Some(a) => astNodeToString(a).map(Some.apply).toValidated
      case None => None.validNel
    }

    implicit val astNodeToCallBodyElement: CheckedAtoB[GenericAstNode, CallBodyElement] = astNodeToAst andThen CheckedAtoB.fromCheck(convertBodyElement _)

    val callBodyValidation: ErrorOr[Option[CallBodyElement]] = ast.getAttributeAsOptional[CallBodyElement]("body").toValidated

    (callableNameValidation, aliasValidation, afterValidation, callBodyValidation) mapN { (name, alias, after, body) =>
      CallElement(name, alias, after, body)
    }
  }

}
