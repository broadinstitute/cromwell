package wdl.draft3.transforms.linking.expression.types

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.TypeEvaluator
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wom.types.{WomArrayType, WomIntegerType, WomType}

object EngineFunctionEvaluators {
  implicit val rangeFunctionEvaluator: TypeEvaluator[Range] = new TypeEvaluator[Range] {
    override def evaluateType(a: Range, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomIntegerType).map(_ => WomArrayType(WomIntegerType))
    }
  }

  private def validateParamType(param: ExpressionElement, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle], expectedType: WomType): ErrorOr[Unit] = {
    param.evaluateType(linkedValues).flatMap { foundType =>
      if (expectedType.isCoerceableFrom(foundType)) { ().validNel } else { s"Invalid parameter '$param'. Expected '${expectedType.toDisplayString}' but got '${foundType.toDisplayString}'".invalidNel }
    }
  }
}
