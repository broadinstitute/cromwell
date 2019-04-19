package wdl.transforms.base.ast2wdlom

import common.Checked
import common.transforms.CheckedAtoB
import common.validation.Checked._
import wdl.model.draft3.elements.ExpressionElement.{StringEscapeSequence, StringLiteral, StringPiece, StringPlaceholder}
import wdl.model.draft3.elements._

object AstNodeToStringPiece {
  def astNodeToStringPiece(expressionProcessor: Option[CheckedAtoB[GenericAstNode, ExpressionElement]]): CheckedAtoB[GenericAstNode, StringPiece] = CheckedAtoB.fromCheck("read string piece") {

    case simple: GenericTerminal if simple.getTerminalStr == "string" => StringLiteral(simple.getSourceString).validNelCheck
    case escape: GenericTerminal if escape.getTerminalStr == "escape" => StringEscapeSequence.parseEscapeSequence(escape.getSourceString).toEither
    case expr: GenericAst if expr.getName == "ExpressionPlaceholder" =>
      expressionProcessor match {
        case Some(processor) => expr.getAttributeAs[ExpressionElement]("expr")(processor) map StringPlaceholder: Checked[StringPlaceholder]
        case None => "String placeholders are not allowed in static strings".invalidNelCheck
      }


    case otherTerminal: GenericTerminal => s"Unexpected parse tree. Expected string piece but found Terminal '${otherTerminal.getTerminalStr}' (${otherTerminal.getSourceString})".invalidNelCheck
    case otherAst: GenericAst => s"Unexpected parse tree. Expected string piece but found AST ${otherAst.getName}".invalidNelCheck
  }
}

object AstNodeToStaticString {
  def astNodeToStaticStringElement(): CheckedAtoB[GenericAstNode, StaticString] = CheckedAtoB.fromCheck("convert AstNode to StaticString") {
    case a: GenericAst if a.getName == "StaticString" =>

      implicit val astNodeToStringPiece = AstNodeToStringPiece.astNodeToStringPiece(None)

      if (a.getAttributes.contains("value")) {
        a.getAttributeAsVector[StringPiece]("value") map { pieces =>
          val unescaped = pieces map {
            case StringLiteral(s) => s
            case e: StringEscapeSequence => e.unescape
            case _: StringPlaceholder => "String placeholders are not allowed in static strings".invalidNelCheck
          }
          StaticString(unescaped.mkString(""))
        }
      } else StaticString("").validNelCheck

    case other => s"Bad value $other".invalidNelCheck
  }
}
