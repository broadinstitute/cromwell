package wdl.draft3.transforms.linking.expression

import cats.syntax.apply._
import cats.syntax.traverse._
import cats.syntax.validated._
import cats.instances.list._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.model.draft3.graph.expression.ValueEvaluator
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wom.expression.IoFunctionSet
import wom.values.{WomArray, WomMap, WomObject, WomPair, WomString, WomValue}

package object values {
  implicit val primitiveValueEvaluator: ValueEvaluator[PrimitiveLiteralExpressionElement] = new ValueEvaluator[PrimitiveLiteralExpressionElement] {
    override def evaluateValue(a: PrimitiveLiteralExpressionElement,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomValue] =
      a.value.validNel
  }

  implicit val stringLiteralEvaluator: ValueEvaluator[StringLiteral] = new ValueEvaluator[StringLiteral] {
    override def evaluateValue(a: StringLiteral,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomValue] =
      WomString(a.value).validNel
  }

  implicit val identifierLookupEvaluator: ValueEvaluator[IdentifierLookup] = new ValueEvaluator[IdentifierLookup] {
    override def evaluateValue(a: IdentifierLookup,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomValue] = {
      inputs.get(a.identifier) match {
        case Some(value) => value.validNel
        case None => s"No suitable input for identifier lookup '${a.identifier}' amongst {${inputs.keys.mkString(", ")}}".invalidNel
      }
    }
  }

  implicit val objectLiteralEvaluator: ValueEvaluator[ObjectLiteral] = new ValueEvaluator[ObjectLiteral] {
    override def evaluateValue(a: ObjectLiteral,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomValue] = {

      val evaluated: ErrorOr[List[(String, WomValue)]] = a.elements.toList traverse { case (key, value) =>
        value.evaluateValue(inputs, ioFunctionSet, linkedValues).map(key -> _)
      }

      evaluated map { l => WomObject(l.toMap) }
    }
  }

  implicit val mapLiteralEvaluator: ValueEvaluator[MapLiteral] = new ValueEvaluator[MapLiteral] {
    override def evaluateValue(a: MapLiteral,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomValue] = {

      val evaluated: ErrorOr[List[(WomValue, WomValue)]] = a.elements.toList traverse { case (key, value) =>
        (key.evaluateValue(inputs, ioFunctionSet, linkedValues),
          value.evaluateValue(inputs, ioFunctionSet, linkedValues)) mapN { (key, value) => key -> value}
      }

      evaluated map { kvps => WomMap(kvps.toMap) }
    }
  }

  implicit val arrayLiteralEvaluator: ValueEvaluator[ArrayLiteral] = new ValueEvaluator[ArrayLiteral] {
    override def evaluateValue(a: ArrayLiteral,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomValue] = {

      val evaluated: ErrorOr[Seq[WomValue]] = a.elements.toList traverse { entry =>
        entry.evaluateValue(inputs, ioFunctionSet, linkedValues)
      }

      evaluated map WomArray.apply
    }
  }

  implicit val pairLiteralEvaluator: ValueEvaluator[PairLiteral] = new ValueEvaluator[PairLiteral] {
    override def evaluateValue(a: PairLiteral,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomValue] = {

      (a.left.evaluateValue(inputs, ioFunctionSet, linkedValues), a.right.evaluateValue(inputs, ioFunctionSet, linkedValues)) mapN { (left, right) =>
        WomPair(left, right)
      }
    }
  }

  implicit val expressionEvaluator: ValueEvaluator[ExpressionElement] = new ValueEvaluator[ExpressionElement] {
    override def evaluateValue(a: ExpressionElement, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomValue] = a match {
      case primitive: PrimitiveLiteralExpressionElement => primitive.evaluateValue(inputs, ioFunctionSet, linkedValues)
      case s: StringLiteral => s.evaluateValue(inputs, ioFunctionSet, linkedValues)
      case id: IdentifierLookup => id.evaluateValue(inputs, ioFunctionSet, linkedValues)
      case o: ObjectLiteral => o.evaluateValue(inputs, ioFunctionSet, linkedValues)
      case m: MapLiteral => m.evaluateValue(inputs, ioFunctionSet, linkedValues)
      case a: ArrayLiteral => a.evaluateValue(inputs, ioFunctionSet, linkedValues)
      case p: PairLiteral => p.evaluateValue(inputs, ioFunctionSet, linkedValues)
      case other => s"Unable to process ${other.getClass.getSimpleName}: No evaluateValue exists for that type.".invalidNel
    }
  }
}
