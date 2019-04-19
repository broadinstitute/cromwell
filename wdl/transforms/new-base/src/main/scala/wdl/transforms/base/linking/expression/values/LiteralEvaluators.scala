package wdl.transforms.base.linking.expression.values

import cats.syntax.apply._
import cats.syntax.traverse._
import cats.syntax.validated._
import cats.instances.list._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.{EvaluatedValue, ForCommandInstantiationOptions, ValueEvaluator}
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wom.expression.IoFunctionSet
import wom.values.{WomArray, WomMap, WomObject, WomPair, WomString, WomValue}


object LiteralEvaluators {

  implicit val primitiveValueEvaluator: ValueEvaluator[PrimitiveLiteralExpressionElement] = new ValueEvaluator[PrimitiveLiteralExpressionElement] {
    override def evaluateValue(a: PrimitiveLiteralExpressionElement,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {
      EvaluatedValue(a.value, Seq.empty).validNel
    }
  }

  implicit val stringLiteralEvaluator: ValueEvaluator[StringLiteral] = new ValueEvaluator[StringLiteral] {
    override def evaluateValue(a: StringLiteral,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomString]] = {
      EvaluatedValue(WomString(a.value), Seq.empty).validNel
    }
  }

  implicit val stringExpressionEvaluator: ValueEvaluator[StringExpression] = new ValueEvaluator[StringExpression] {
    override def evaluateValue(a: StringExpression,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {
      val evaluatedPieces = a.pieces.toList.traverse {
        case e: StringPlaceholder => e.expr.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case s: StringLiteral => s.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case e: StringEscapeSequence => EvaluatedValue(WomString(e.unescape), Seq.empty).validNel
      }

      evaluatedPieces map { pieces =>
        EvaluatedValue(WomString(pieces.map(_.value.valueString).mkString("")), pieces.flatMap(_.sideEffectFiles))
      }
    }
  }

  implicit val objectLiteralEvaluator: ValueEvaluator[ObjectLiteral] = new ValueEvaluator[ObjectLiteral] {
    override def evaluateValue(a: ObjectLiteral,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {

      val evaluated: ErrorOr[List[(String, EvaluatedValue[_])]] = a.elements.toList traverse { case (key, value) =>
        value.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions).map(key -> _)
      }

      evaluated map { mapping =>
        val value = mapping.map(entry => entry._1 -> entry._2.value).toMap
        val sideEffectFiles = mapping.flatMap(entry => entry._2.sideEffectFiles)
        EvaluatedValue(WomObject(value), sideEffectFiles)
      }
    }
  }

  implicit val mapLiteralEvaluator: ValueEvaluator[MapLiteral] = new ValueEvaluator[MapLiteral] {
    override def evaluateValue(a: MapLiteral,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {

      val evaluated: ErrorOr[List[(EvaluatedValue[_], EvaluatedValue[_])]] = a.elements.toList traverse { case (key, value) =>
        (key.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions),
          value.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) mapN { (key, value) => key -> value}
      }

      evaluated map { kvps =>
        val value = kvps.map(entry => entry._1.value -> entry._2.value).toMap
        val sideEffectFiles = kvps.flatMap(entry => entry._1.sideEffectFiles ++ entry._2.sideEffectFiles)
        EvaluatedValue(WomMap(value.toMap), sideEffectFiles)
      }
    }
  }

  implicit val arrayLiteralEvaluator: ValueEvaluator[ArrayLiteral] = new ValueEvaluator[ArrayLiteral] {
    override def evaluateValue(a: ArrayLiteral,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {

      val evaluated: ErrorOr[Seq[EvaluatedValue[_]]] = a.elements.toList traverse { entry =>
        entry.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
      }

      evaluated map { evaluateds =>
        val value = evaluateds.map(_.value)
        val sideEffectFiles = evaluateds.flatMap(_.sideEffectFiles)
        EvaluatedValue(WomArray(value), sideEffectFiles)
      }
    }
  }

  implicit val pairLiteralEvaluator: ValueEvaluator[PairLiteral] = new ValueEvaluator[PairLiteral] {
    override def evaluateValue(a: PairLiteral,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {

      (a.left.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions),
        a.right.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) mapN { (left, right) =>

        EvaluatedValue(WomPair(left.value, right.value), left.sideEffectFiles ++ right.sideEffectFiles)
      }
    }
  }
}
