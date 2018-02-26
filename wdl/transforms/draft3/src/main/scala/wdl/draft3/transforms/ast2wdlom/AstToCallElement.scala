package wdl.draft3.transforms.ast2wdlom

import cats.syntax.validated._
import cats.syntax.either._
import cats.syntax.apply._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.Ast
import wdl.model.draft3.elements.CallElement
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wdl.model.draft3.elements.ExpressionElement.KvPair

object AstToCallElement {
  def convert(ast: Ast): ErrorOr[CallElement] = {
    val callableNameValidation: ErrorOr[String] = astNodeToString(ast.getAttribute("task")).toValidated
    val aliasValidation: ErrorOr[Option[String]] = Option(ast.getAttribute("alias")) match {
      case Some(a) => astNodeToString(a).map(Some.apply).toValidated
      case None => None.validNel
    }
    val callBodyValidation: ErrorOr[Vector[KvPair]] = Option(ast.getAttributeAsVector[KvPair]("body")) match {
      case Some(p) => p.toValidated
      case None => Vector.empty.validNel
    }

    (callableNameValidation, aliasValidation, callBodyValidation) mapN { (name, alias, body) =>
      CallElement(name, alias, body)
    }
  }

}
