package wdl.draft3.transforms.ast2wdlom

import cats.syntax.validated._
import cats.syntax.either._
import cats.instances.either._
import cats.syntax.apply._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.{Ast, AstNode}
import wdl.model.draft3.elements.{CallBodyElement, CallElement}
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wdl.model.draft3.elements.ExpressionElement.KvPair

object AstToCallElement {
  def convert(ast: Ast): ErrorOr[CallElement] = {

    val callableNameValidation: ErrorOr[String] = astNodeToString(ast.getAttribute("task")).toValidated

    val aliasValidation: ErrorOr[Option[String]] = Option(ast.getAttribute("alias")) match {
      case Some(a) => astNodeToString(a).map(Some.apply).toValidated
      case None => None.validNel
    }

    implicit val astNodeToCallBodyElement: CheckedAtoB[AstNode, CallBodyElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToCallBodyElement.convert)

    val callBodyValidation: ErrorOr[Option[CallBodyElement]] = ast.getAttributeAsOptional[CallBodyElement]("body").toValidated

    (callableNameValidation, aliasValidation, callBodyValidation) mapN { (name, alias, body) =>
      CallElement(name, alias, body)
    }
  }

  object AstToCallBodyElement {
    def convert(a: Ast): Checked[CallBodyElement] = {
      a.getAttributeAsVector[KvPair]("inputs") map CallBodyElement
    }
  }
}
