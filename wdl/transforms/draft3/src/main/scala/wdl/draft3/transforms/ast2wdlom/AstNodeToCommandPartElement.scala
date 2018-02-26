package wdl.draft3.transforms.ast2wdlom

import cats.syntax.validated._
import cats.syntax.either._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser._
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wdl.model.draft3.elements.{CommandPartElement, ExpressionElement}
import wdl.model.draft3.elements.CommandPartElement.{PlaceholderCommandPartElement, StringCommandPartElement}

object AstNodeToCommandPartElement {
  def convert(a: AstNode): ErrorOr[CommandPartElement] = a match {
    case t: Terminal => astNodeToString(t).toValidated map StringCommandPartElement
    case a: Ast => a.getAttributeAs[ExpressionElement]("expr").toValidated map PlaceholderCommandPartElement
    case other => s"Conversion for $other not supported".invalidNel
  }
}
