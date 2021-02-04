package wdl.model.draft3.graph.expression

import common.validation.ErrorOr.ErrorOr
import simulacrum.typeclass
import wdl.model.draft3.elements.ExpressionElement
import wom.CommandSetupSideEffectFile
import wom.expression.IoFunctionSet
import wom.values.WomValue

@typeclass
trait ValueEvaluator[A] {
  /**
    * Evaluate a value from an A
    * @param a The A to evaluate
    * @param inputs Evaluation inputs
    * @param ioFunctionSet IO functions to use
    * @param forCommandInstantiationOptions Supplied only if we're evaluating this A as part of command instantiation.
    * @return An evaluated value set - the value itself and any files which were produced as part of the evaluation.
    */
  def evaluateValue(a: A,
                    inputs: Map[String, WomValue],
                    ioFunctionSet: IoFunctionSet,
                    forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                   (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]]
}

final case class ForCommandInstantiationOptions(valueMapper: WomValue => WomValue)
final case class EvaluatedValue[A <: WomValue](value: A with WomValue, sideEffectFiles: Seq[CommandSetupSideEffectFile])
