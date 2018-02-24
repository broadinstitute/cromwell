package wdl.draft3.transforms.linking.expression

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement.{IdentifierLookup, PrimitiveLiteralExpressionElement}
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook, UnlinkedIdentifierHook}
import wdl.model.draft3.graph.expression.ValueEvaluator
import wom.values.WomValue

package object values {
  implicit object primitiveValueEvaluator extends ValueEvaluator[PrimitiveLiteralExpressionElement] {
    override def evaluateValue(a: PrimitiveLiteralExpressionElement,
                               inputs: Map[String, WomValue],
                               linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomValue] =
      a.value.validNel
  }

  implicit object identifierLookupEvaluator extends ValueEvaluator[IdentifierLookup] {
    override def evaluateValue(a: IdentifierLookup,
                               inputs: Map[String, WomValue],
                               linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomValue] = {
      linkedValues.collectFirst {
        case (UnlinkedIdentifierHook(id), generatedValueHandle) if a.identifier == id => generatedValueHandle.linkableName
      } flatMap { inputs.get } match {
        case Some(value) => value.validNel
        case None => s"No suitable input for identifier lookup '${a.identifier}' amongst {${inputs.keys.mkString(", ")}}".invalidNel
      }
    }
  }

}
