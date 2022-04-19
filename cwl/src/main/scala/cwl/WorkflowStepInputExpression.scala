package cwl

import cats.syntax.validated._
import wom.expression.{FileEvaluation, IoFunctionSet}
import wom.types._
import wom.values._

final case class WorkflowStepInputExpression(inputName: String,
                                             valueFrom: StringOrExpression,
                                             override val cwlExpressionType: WomType,
                                             inputs: Set[String],
                                             override val expressionLib: ExpressionLib) extends CwlWomExpression {

  override def sourceString = inputName

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet) = {
    valueFrom match {
      // If valueFrom is a constant string value, use this as the value for this input parameter.
      // TODO: need to handle case where this is a parameter reference, it currently looks like a String to us!
      case StringOrExpression.String(value) => WomString(value).validNel

      /*
       * If valueFrom is a parameter reference or expression, it must be evaluated to yield the actual value to be assiged to the input field.
       *
       * The self value of in the parameter reference or expression must be the value of the parameter(s) specified in the source field,
       * or null if there is no source field.
       *
       * The value of inputs in the parameter reference or expression must be the input object to the workflow step after
       * assigning the source values, applying default, and then scattering. The order of evaluating valueFrom among step
       * input parameters is undefined and the result of evaluating valueFrom on a parameter must not be visible to
       * evaluation of valueFrom on other parameters.
       */
      case StringOrExpression.Expression(expression) =>
        //used to determine the value of "self" as expected by CWL Spec
        def selfValue = inputValues.get(inputName) match {
          case Some(value) => value
          case None => WomOptionalValue(WomNothingType, None)
        }

        val parameterContext = ParameterContext(ioFunctionSet, expressionLib, inputValues, selfValue)

        expression.fold(EvaluateExpression).apply(parameterContext)
      case oh => throw new Exception(s"Programmer Error! Unexpected case match: $oh")
    }
  }

  //this is input, so not producing any output files
  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType) =
    Set.empty[FileEvaluation].validNel
}

