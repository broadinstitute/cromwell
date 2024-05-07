package wdl.transforms.base.ast2wdlom

import cats.data.NonEmptyList
import cats.syntax.validated._
import cats.syntax.either._
import cats.syntax.apply._
import cats.syntax.traverse._
import cats.instances.vector._
import cats.instances.either._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.Checked._
import common.validation.ErrorOr._
import wdl.model.draft3.elements.ExpressionElement.{PrimitiveLiteralExpressionElement, StringExpression, StringLiteral}
import wdl.model.draft3.elements._

object AstNodeToPlaceholderAttributeSet {

  def attributeKvpConverter(implicit
    astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement]
  ): CheckedAtoB[GenericAstList, PlaceholderAttributeSet] = {
    val singleElement = astNodeToAst andThen CheckedAtoB.fromErrorOr(convertAttributeKvp _)

    def convertAll(as: GenericAstList): ErrorOr[Vector[PlaceholderAttributeElement]] =
      as.astNodeList.toVector.traverse[ErrorOr, PlaceholderAttributeElement](singleElement.run(_).toValidated)

    CheckedAtoB.fromErrorOr(convertAll _) andThen CheckedAtoB.fromCheck(validAttributeSet _)
  }

  private def validAttributeSet(a: Vector[PlaceholderAttributeElement]): Checked[PlaceholderAttributeSet] = {
    def foldFunction(acc: (PlaceholderAttributeSet, List[String]),
                     next: PlaceholderAttributeElement
    ): (PlaceholderAttributeSet, List[String]) =
      (acc, next) match {
        case ((PlaceholderAttributeSet(None, t, f, s), errors), DefaultAttributeElement(d)) =>
          (PlaceholderAttributeSet(Option(d), t, f, s), errors)
        case ((pas @ PlaceholderAttributeSet(Some(d1), _, _, _), errors), DefaultAttributeElement(d2)) =>
          (pas, errors :+ s"""Cannot supply 'default="$d2"' because it is already supplied as 'default="$d1"'""")

        case ((PlaceholderAttributeSet(d, None, f, s), errors), TrueAttributeElement(t)) =>
          (PlaceholderAttributeSet(d, Option(t), f, s), errors)
        case ((pas @ PlaceholderAttributeSet(_, Some(t1), _, _), errors), TrueAttributeElement(t2)) =>
          (pas, errors :+ s"""Cannot supply 'true="$t2"' because it is already supplied as 'true="$t1"'""")

        case ((PlaceholderAttributeSet(d, t, None, s), errors), FalseAttributeElement(f)) =>
          (PlaceholderAttributeSet(d, t, Option(f), s), errors)
        case ((pas @ PlaceholderAttributeSet(_, _, Some(f1), _), errors), FalseAttributeElement(f2)) =>
          (pas, errors :+ s"""Cannot supply 'false="$f2"' because it is already supplied as 'false="$f1"'""")

        case ((PlaceholderAttributeSet(d, t, f, None), errors), SepAttributeElement(s)) =>
          (PlaceholderAttributeSet(d, t, f, Option(s)), errors)
        case ((pas @ PlaceholderAttributeSet(_, _, _, Some(s1)), errors), SepAttributeElement(s2)) =>
          (pas, errors :+ s"""Cannot supply 'sep="$s2"' because it is already supplied as 'sep="$s1"'""")
      }

    val folded = a.foldLeft((PlaceholderAttributeSet.empty, List.empty[String]))(foldFunction)

    folded match {
      case (_, nel) if nel.nonEmpty => Left(NonEmptyList.fromListUnsafe(nel))
      case (PlaceholderAttributeSet(_, Some(t), None, _), _) =>
        s"""Cannot specify 'true="$t"' without also having a 'false="..."' attribute""".invalidNelCheck
      case (PlaceholderAttributeSet(_, None, Some(f), _), _) =>
        s"""Cannot specify 'false="$f"' without also having a 'true="..."' attribute""".invalidNelCheck
      case (valid, _) => valid.validNelCheck
    }
  }

  private def convertAttributeKvp(a: GenericAst)(implicit
    astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement]
  ): ErrorOr[PlaceholderAttributeElement] =
    (placeholderAttributeConstructor(a), placeholderAttributeValue(a)) mapN { (constructor, value) =>
      constructor.apply(value)
    }

  private def placeholderAttributeConstructor(kvpAst: GenericAst): ErrorOr[String => PlaceholderAttributeElement] =
    kvpAst.getAttributeAs[String]("key").toValidated.flatMap {
      case "sep" => (SepAttributeElement.apply _).validNel
      case "default" => (DefaultAttributeElement.apply _).validNel
      case "true" => (TrueAttributeElement.apply _).validNel
      case "false" => (FalseAttributeElement.apply _).validNel
      case other => s"Invalid placeholder attribute: $other".invalidNel
    }

  private def placeholderAttributeValue(
    kvpAst: GenericAst
  )(implicit astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement]): ErrorOr[String] =
    kvpAst.getAttributeAs[ExpressionElement]("value").toValidated.flatMap {
      case StringLiteral(literalValue) => literalValue.validNel
      case StringExpression(pieces) if pieces.length == 1 =>
        pieces.head match {
          case StringLiteral(literalValue) => literalValue.validNel
          case other => s"Cannot use $other as a placeholder attribute. Must be a primitive literal".invalidNel
        }
      case PrimitiveLiteralExpressionElement(primitive) => primitive.valueString.validNel
      case other => s"Cannot use $other as a placeholder attribute. Must be a primitive literal".invalidNel
    }
}
