package wdl.transforms.base.ast2wdlom

import common.validation.Validation._
import cats.syntax.validated._
import cats.syntax.apply._
import cats.syntax.either._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.StaticString
import wom.callable.MetaValueElement._
import wom.callable.{MetaKvPair, MetaValueElement}

import scala.util.{Failure, Try}

object AstNodeToMetaKvPair {
  def astNodeToMetaKvPair(implicit astNodeToStaticString: CheckedAtoB[GenericAstNode, StaticString]): CheckedAtoB[GenericAstNode, MetaKvPair] = {
    CheckedAtoB.fromErrorOr("convert AstNode to MetaKvPair")(convertMetaKvPair)
  }

  private def convertMetaKvPair(astNode: GenericAstNode)
                               (implicit astNodeToStaticString: CheckedAtoB[GenericAstNode, StaticString]): ErrorOr[MetaKvPair] = astNode match {
    case a: GenericAst if a.getName == "MetaKvPair" =>
      val keyValidation: ErrorOr[String] = a.getAttributeAs[String]("key").toValidated
      val valueValidation: ErrorOr[MetaValueElement] = convertMetaValue(a.getAttribute("value"))

      (keyValidation, valueValidation) mapN { (key, value) => MetaKvPair(key, value) }
    case other => s"Expected Ast of type 'MetaKvPair' but got $other".invalidNel
  }

  private def convertMetaValue(astNode: GenericAstNode)
                              (implicit astNodeToStaticString: CheckedAtoB[GenericAstNode, StaticString]): ErrorOr[MetaValueElement] = {
    implicit val recursiveKvPairConversion = CheckedAtoB.fromErrorOr(convertMetaKvPair _)
    astNode match {
      // This is a primitive type, one of {null, boolean, float, int, string}.
      case t: GenericTerminal =>
        (t.getTerminalStr, t.getSourceString) match {
          case ("integer", i) => Try(MetaValueElementInteger(i.toInt)).toErrorOr
          case ("float", f) => Try(MetaValueElementFloat(f.toDouble)).toErrorOr
          case ("boolean", b) => Try(MetaValueElementBoolean(b.toBoolean)).toErrorOr
          case ("string", s) => Try(MetaValueElementString(s)).toErrorOr
          case ("null", _) => MetaValueElementNull.validNel
          case (name, other) => s"No conversion defined for Ast ($name, $other) to MetaValueElement".invalidNel
        }

      case a: GenericAst if a.getName == "StaticString" => astNodeToStaticString(a).toValidated map { staticString => MetaValueElementString.apply(staticString.value) }

      case a: GenericAst if a.getName == "MetaArray" =>
        a.getAttributeAsVectorF[MetaValueElement]("values")(convertMetaValue(_).toEither).toValidated.map(MetaValueElementArray)

      case a: GenericAst if a.getName == "MetaObject" =>
        (for {
          mapKvs <- a.getAttributeAsVector[MetaKvPair]("map")
          asMap = mapKvs.map(kv => kv.key -> kv.value).toMap
        } yield MetaValueElementObject(asMap)).toValidated

      case other =>
        Failure(new Exception(s"No conversion defined for AstNode $other to MetaValueElement")).toErrorOr
    }
  }
}
