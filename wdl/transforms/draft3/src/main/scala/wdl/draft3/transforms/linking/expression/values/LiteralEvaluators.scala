package wdl.draft3.transforms.linking.expression.values

import cats.syntax.apply._
import cats.syntax.traverse._
import cats.syntax.validated._
import cats.instances.list._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.model.draft3.graph.expression.ValueEvaluator
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wom.expression.IoFunctionSet
import wom.values.{WomArray, WomMap, WomObject, WomPair, WomString, WomValue}


object LiteralEvaluators {
  implicit val primitiveValueEvaluator: ValueEvaluator[PrimitiveLiteralExpressionElement] = new ValueEvaluator[PrimitiveLiteralExpressionElement] {
    override def evaluateValue(a: PrimitiveLiteralExpressionElement,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] =
      a.value.validNel
  }

  implicit val stringLiteralEvaluator: ValueEvaluator[StringLiteral] = new ValueEvaluator[StringLiteral] {
    override def evaluateValue(a: StringLiteral,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] =
      WomString(a.value).validNel
  }

  implicit val objectLiteralEvaluator: ValueEvaluator[ObjectLiteral] = new ValueEvaluator[ObjectLiteral] {
    override def evaluateValue(a: ObjectLiteral,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {

      val evaluated: ErrorOr[List[(String, WomValue)]] = a.elements.toList traverse { case (key, value) =>
        value.evaluateValue(inputs, ioFunctionSet).map(key -> _)
      }

      evaluated map { l => WomObject(l.toMap) }
    }
  }

  implicit val mapLiteralEvaluator: ValueEvaluator[MapLiteral] = new ValueEvaluator[MapLiteral] {
    override def evaluateValue(a: MapLiteral,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {

      val evaluated: ErrorOr[List[(WomValue, WomValue)]] = a.elements.toList traverse { case (key, value) =>
        (key.evaluateValue(inputs, ioFunctionSet),
          value.evaluateValue(inputs, ioFunctionSet)) mapN { (key, value) => key -> value}
      }

      evaluated map { kvps => WomMap(kvps.toMap) }
    }
  }

  implicit val arrayLiteralEvaluator: ValueEvaluator[ArrayLiteral] = new ValueEvaluator[ArrayLiteral] {
    override def evaluateValue(a: ArrayLiteral,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {

      val evaluated: ErrorOr[Seq[WomValue]] = a.elements.toList traverse { entry =>
        entry.evaluateValue(inputs, ioFunctionSet)
      }

      evaluated map WomArray.apply
    }
  }

  implicit val pairLiteralEvaluator: ValueEvaluator[PairLiteral] = new ValueEvaluator[PairLiteral] {
    override def evaluateValue(a: PairLiteral,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {

      (a.left.evaluateValue(inputs, ioFunctionSet),
        a.right.evaluateValue(inputs, ioFunctionSet)) mapN { (left, right) => WomPair(left, right) }
    }
  }
}
