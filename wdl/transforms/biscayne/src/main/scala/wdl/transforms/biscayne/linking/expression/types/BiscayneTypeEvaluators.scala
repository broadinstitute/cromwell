package wdl.transforms.biscayne.linking.expression.types

import cats.syntax.validated._
import common.validation.ErrorOr._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.{AsMap, AsPairs, CollectByKey, Keys}
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.model.draft3.graph.expression.TypeEvaluator
import wdl.transforms.base.linking.expression.types.EngineFunctionEvaluators.validateParamType
import wom.types._

object BiscayneTypeEvaluators {
  implicit val keysFunctionEvaluator: TypeEvaluator[Keys] = new TypeEvaluator[Keys] {
    override def evaluateType(a: Keys, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomMapType(WomAnyType, WomAnyType)) flatMap {
        case WomMapType(keyType, _) => WomArrayType(keyType).validNel
        case other => s"Cannot invoke 'keys' on type '${other.stableName}'. Expected a map".invalidNel
      }
    }
  }

  implicit val asMapFunctionEvaluator: TypeEvaluator[AsMap] = new TypeEvaluator[AsMap] {
    override def evaluateType(a: AsMap, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomArrayType(WomPairType(WomAnyType, WomAnyType))) flatMap {
        case WomArrayType(WomPairType(x: WomPrimitiveType, y)) => WomMapType(x, y).validNel
        case other @ WomArrayType(WomPairType(x, _)) => s"Cannot invoke 'as_map' on type ${other.stableName}. Map keys must be primitive but got '${x.stableName}'".invalidNel
        case other => s"Cannot invoke 'as_map' on type '${other.stableName}'. Expected an array of pairs".invalidNel
      }
    }
  }

  implicit val asPairsFunctionEvaluator: TypeEvaluator[AsPairs] = new TypeEvaluator[AsPairs] {
    override def evaluateType(a: AsPairs, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomMapType(WomAnyType, WomAnyType)) flatMap {
        case WomMapType(x, y) => WomArrayType(WomPairType(x, y)).validNel
        case other => s"Cannot invoke 'as_pairs' on type '${other.stableName}'. Expected a map".invalidNel
      }
    }
  }

  implicit val collectByKeyFunctionEvaluator: TypeEvaluator[CollectByKey] = new TypeEvaluator[CollectByKey] {
    override def evaluateType(a: CollectByKey, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomArrayType(WomPairType(WomAnyType, WomAnyType))) flatMap {
        case WomArrayType(WomPairType(x: WomPrimitiveType, y)) => WomMapType(x, WomArrayType(y)).validNel
        case other@WomArrayType(WomPairType(x, _)) => s"Cannot invoke 'collect_by_key' on type ${other.stableName}. Map keys must be primitive but got '${x.stableName}'".invalidNel
        case other => s"Cannot invoke 'collect_by_key' on type '${other.stableName}'. Expected an array of pairs".invalidNel
      }
    }
  }
}
