package wdl.transforms.biscayne.linking.expression.values

import cats.data.NonEmptyList
import cats.syntax.validated._
import cats.syntax.traverse._
import cats.instances.list._
import common.validation.ErrorOr._
import common.collections.EnhancedCollections._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.{EvaluatedValue, ForCommandInstantiationOptions, ValueEvaluator}
import wdl.transforms.base.linking.expression.values.EngineFunctionEvaluators.processValidatedSingleValue
import wom.expression.IoFunctionSet
import wom.types._
import wom.values.{WomArray, WomMap, WomOptionalValue, WomPair, WomValue}
import wom.types.coercion.defaults._

object BiscayneValueEvaluators {

  implicit val noneLiteralEvaluator: ValueEvaluator[NoneLiteralElement.type] = new ValueEvaluator[ExpressionElement.NoneLiteralElement.type] {
    override def evaluateValue(a: ExpressionElement.NoneLiteralElement.type, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])(implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {
      EvaluatedValue(
        value = WomOptionalValue(WomNothingType, None),
        sideEffectFiles = Seq.empty).validNel
    }
  }


  implicit val asMapFunctionEvaluator: ValueEvaluator[AsMap] = new ValueEvaluator[AsMap] {
    override def evaluateValue(a: AsMap, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {
      processValidatedSingleValue[WomArray, WomMap](expressionValueEvaluator.evaluateValue(a.param, inputs, ioFunctionSet, forCommandInstantiationOptions)(expressionValueEvaluator)) {
        case WomArray(WomArrayType(WomPairType(_: WomPrimitiveType, _)), values) =>
          val validPairs: ErrorOr[List[(WomValue, WomValue)]] = values.toList traverse {
            case WomPair(l, r) => (l, r).validNel
            case other => s"Unexpected array element. Expected a Pair[X, Y] but array contained ${other.toWomString}]".invalidNel
          } leftMap {
            errors => NonEmptyList.fromListUnsafe(errors.toList.distinct)
          }
          validPairs flatMap { pairs =>
            val grouped = pairs.groupBy(_._1)
            val tooManyKeyErrors = grouped collect {
              case (name, list) if list.length != 1 => s"keys can only appear once but ${name.toWomString} appeared ${list.size} times."
            }
            if (tooManyKeyErrors.isEmpty) {
              val pairs = grouped map { case (key, value) => (key -> value.head._2) }

              EvaluatedValue(WomMap(pairs), Seq.empty).validNel
            } else {
              s"Cannot evaluate 'as_map' with duplicated keys: ${tooManyKeyErrors.mkString(", ")}".invalidNel
            }
          }

        case WomArray(womType@WomArrayType(WomPairType(x, _)), _) => s"Cannot evaluate 'as_map' on type ${womType.stableName}. Keys must be primitive but got ${x.stableName}.".invalidNel
        case other => s"Invalid call of 'as_map' on parameter of type '${other.womType.stableName}' (expected Array[Pair[X, Y]])".invalidNel
      } (coercer = WomArrayType(WomPairType(WomAnyType, WomAnyType)))
    }
  }

  implicit val keysFunctionEvaluator: ValueEvaluator[Keys] = new ValueEvaluator[Keys] {
    override def evaluateValue(a: Keys,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomArray]] = {

      processValidatedSingleValue[WomMap, WomArray](expressionValueEvaluator.evaluateValue(a.param, inputs, ioFunctionSet, forCommandInstantiationOptions)(expressionValueEvaluator)) {
        case WomMap(WomMapType(keyType, _), values) => EvaluatedValue(WomArray(WomArrayType(keyType), values.keys.toList), Seq.empty).validNel
        case other => s"Invalid call of 'keys' on parameter of type '${other.womType.stableName}' (expected Map[X, Y])".invalidNel
      }
    }
  }

  implicit val asPairsFunctionEvaluator: ValueEvaluator[AsPairs] = new ValueEvaluator[AsPairs] {
    override def evaluateValue(a: AsPairs, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {
      processValidatedSingleValue[WomMap, WomArray](expressionValueEvaluator.evaluateValue(a.param, inputs, ioFunctionSet, forCommandInstantiationOptions)(expressionValueEvaluator)) {
        case WomMap(WomMapType(keyType, valueType), values) =>
          val validPairs: List[WomPair] = values.toList map {
            case (l, r) => WomPair(l, r)
          }
          EvaluatedValue(WomArray(WomArrayType(WomPairType(keyType, valueType)), validPairs), Seq.empty).validNel

        case other => s"Invalid call of 'as_pairs' on parameter of type '${other.womType.stableName}' (expected Map[X, Y])".invalidNel
      }
    }
  }

  implicit val collectByKeyFunctionEvaluator: ValueEvaluator[CollectByKey] = new ValueEvaluator[CollectByKey] {
    override def evaluateValue(a: CollectByKey, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {
      processValidatedSingleValue[WomArray, WomMap](expressionValueEvaluator.evaluateValue(a.param, inputs, ioFunctionSet, forCommandInstantiationOptions)(expressionValueEvaluator)) {
        case WomArray(WomArrayType(WomPairType(_: WomPrimitiveType, _)), values) =>
          val validPairs: ErrorOr[List[(WomValue, WomValue)]] = values.toList traverse {
            case WomPair(l, r) => (l, r).validNel
            case other => s"Unexpected array element. Expected a Pair[X, Y] but array contained ${other.toWomString}]".invalidNel
          }
          validPairs flatMap { kvpairs =>
            val grouped: Map[WomValue, WomArray] = kvpairs.groupBy(_._1).safeMapValues(v => WomArray(v.map(_._2)))
            EvaluatedValue(WomMap(grouped), Seq.empty).validNel

          }

        case WomArray(womType@WomArrayType(WomPairType(x, _)), _) => s"Cannot evaluate 'collect_by_key' on type ${womType.stableName}. Keys must be primitive but got ${x.stableName}.".invalidNel
        case other => s"Invalid call of 'collect_by_key' on parameter of type '${other.womType.stableName}' (expected Map[X, Y])".invalidNel
      } (coercer = WomArrayType(WomPairType(WomAnyType, WomAnyType)))
    }
  }
}
