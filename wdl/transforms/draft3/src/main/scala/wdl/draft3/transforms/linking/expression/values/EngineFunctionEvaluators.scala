package wdl.draft3.transforms.linking.expression.values

import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.ValueEvaluator
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wom.expression.IoFunctionSet
import wom.types.{WomArrayType, WomIntegerType}
import wom.values.{WomArray, WomInteger, WomValue}

object EngineFunctionEvaluators {
  implicit val rangeFunctionEvaluator: ValueEvaluator[Range] = new ValueEvaluator[Range] {
    override def evaluateValue(a: Range, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      a.param.evaluateValue(inputs, ioFunctionSet) flatMap {
        case WomInteger(i) => WomArray(WomArrayType(WomIntegerType, guaranteedNonEmpty = i > 0), (0 until i).map(WomInteger)).validNel
        case other => s"Could not evaluate $a: Expected integer argument but got ${other.womType}".invalidNel
      }
    }
  }
}
