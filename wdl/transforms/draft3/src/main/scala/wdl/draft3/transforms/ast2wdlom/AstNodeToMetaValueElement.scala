package wdl.draft3.transforms.ast2wdlom

import cats.syntax.either._
import cats.syntax.validated._
import common.validation.Validation._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.{Ast, AstNode, Terminal}
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wom.callable.{MetaKvPair, MetaValueElement}
import wom.callable.MetaValueElement._

import scala.util.{Failure, Try}

object AstNodeToMetaValueElement {

  def convert(astNode: AstNode): ErrorOr[MetaValueElement] = astNode match {

      // This is a primitive type, one of {null, boolean, float, int, string}.
    case t: Terminal =>
      (t.getTerminalStr, t.getSourceString) match {
        case ("integer", i) => Try(MetaValueElementInteger(i.toInt)).toErrorOr
        case ("float", f) => Try(MetaValueElementFloat(f.toDouble)).toErrorOr
        case ("boolean", b) => Try(MetaValueElementBoolean(b.toBoolean)).toErrorOr
        case ("string", s) => MetaValueElementString(s).validNel
        case ("null", _) => MetaValueElementNull.validNel
        case (name,other) => s"No conversion defined for Ast ($name, $other) to MetaValueElement".invalidNel
      }

    case a: Ast if a.getName == "MetaArray" =>
      a.getAttributeAsVector[MetaValueElement]("values").toValidated.map(MetaValueElementArray)

    case a: Ast if a.getName == "MetaObject" =>
      (for {
        mapKvs <- a.getAttributeAsVector[MetaKvPair]("map")
        asMap = mapKvs.map(kv => kv.key -> kv.value).toMap
       } yield MetaValueElementObject(asMap)).toValidated

    case other =>
      Failure(new Exception(s"No conversion defined for AstNode $other to MetaValueElement")).toErrorOr
  }
}
