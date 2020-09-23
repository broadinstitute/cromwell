package wdl.transforms.base.linking.expression.types

import cats.syntax.apply._
import cats.syntax.validated._
import cats.syntax.traverse._
import cats.instances.list._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.model.draft3.graph.expression.TypeEvaluator
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wom.types._

object LiteralEvaluators {
  implicit val primitiveTypeEvaluator: TypeEvaluator[PrimitiveLiteralExpressionElement] = new TypeEvaluator[PrimitiveLiteralExpressionElement] {
    override def evaluateType(a: PrimitiveLiteralExpressionElement, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] =
      a.value.womType.validNel
  }

  implicit val noneLiteralTypeEvaluator: TypeEvaluator[NoneLiteralElement.type] = new TypeEvaluator[NoneLiteralElement.type] {
    override def evaluateType(a: NoneLiteralElement.type, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] =
      WomOptionalType(WomNothingType).validNel
  }

  implicit val objectLiteralTypeEvaluator: TypeEvaluator[ObjectLiteral] = new TypeEvaluator[ObjectLiteral] {
    override def evaluateType(a: ObjectLiteral, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] = WomObjectType.validNel
  }

  implicit val stringLiteralTypeEvaluator: TypeEvaluator[StringLiteral] = new TypeEvaluator[StringLiteral] {
    override def evaluateType(a: StringLiteral, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] = WomStringType.validNel
  }

  implicit val stringExpressionTypeEvaluator: TypeEvaluator[StringExpression] = new TypeEvaluator[StringExpression] {
    override def evaluateType(a: StringExpression, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] = WomStringType.validNel
  }

  implicit val mapLiteralTypeEvaluator: TypeEvaluator[MapLiteral] = new TypeEvaluator[MapLiteral] {
    override def evaluateType(a: MapLiteral, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] = {

      val keyTypes = a.elements.keySet.toList.traverse { x: ExpressionElement => x.evaluateType(linkedValues) }
      val commonKeyType: ErrorOr[WomType] = keyTypes.map(WomType.homogeneousTypeFromTypes)

      val valueTypes = a.elements.values.toList.traverse { y: ExpressionElement => y.evaluateType(linkedValues) }
      val commonValueType: ErrorOr[WomType] = valueTypes.map(WomType.homogeneousTypeFromTypes)

      (commonKeyType, commonValueType) mapN { (k, v) => WomMapType(k, v) }
    }
  }

  implicit val arrayLiteralTypeEvaluator: TypeEvaluator[ArrayLiteral] = new TypeEvaluator[ArrayLiteral] {
    override def evaluateType(a: ArrayLiteral, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] = {

      val types = a.elements.toList.traverse { x: ExpressionElement => x.evaluateType(linkedValues) }
      val commonType: ErrorOr[WomType] = types.map(WomType.homogeneousTypeFromTypes)

      commonType.map(WomArrayType.apply(_, guaranteedNonEmpty = a.elements.nonEmpty))
    }
  }

  implicit val pairLiteralTypeEvaluator: TypeEvaluator[PairLiteral] = new TypeEvaluator[PairLiteral] {
    override def evaluateType(a: PairLiteral, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] = {

      val leftType = a.left.evaluateType(linkedValues)
      val rightType = a.right.evaluateType(linkedValues)

      (leftType, rightType) mapN { (l, r) => WomPairType(l, r) }
    }
  }
}
