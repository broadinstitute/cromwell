package cwl

import cats.data.NonEmptyList
import cats.syntax.option._
import cats.syntax.validated._
import cats.syntax.either._
import cats.syntax.traverse._
import cats.instances.list._
import common.Checked
import common.validation.ErrorOr.ErrorOr
import cwl.WorkflowStepInput.InputSource
import cwl.command.ParentName
import wom.expression.IoFunctionSet
import wom.types._
import wom.values._
import mouse.all._

final case class WorkflowStepInputExpression(input: WorkflowStepInput,
                                             override val cwlExpressionType: WomType,
                                             graphInputs: Set[String],
                                             override val expressionLib: ExpressionLib)(implicit parentName: ParentName) extends CwlWomExpression {

  override def sourceString = input.toString

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet) = {
    def lookupValue(key: String): Checked[WomValue] =
      inputValues.
        get(key).
        toRight(s"source value $key not found in input values ${inputValues.mkString(", ")}" |> NonEmptyList.one)

    (input.valueFrom, input.source) match {
      case (None, Some(WorkflowStepInputSource.String(id))) =>
        inputValues.
          get(FullyQualifiedName(id).id).
          toValidNel(s"could not find id $id in typeMap\n${inputValues.mkString("\n")}\nwhen evaluating $input.  Graph Inputs were ${graphInputs.mkString("\n")}")

        // If valueFrom is a constant string value, use this as the value for this input parameter.
        // TODO: need to handle case where this is a parameter reference
      case (Some(StringOrExpression.String(value)), None) => WomString(value).validNel

      /**
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
      case (Some(StringOrExpression.Expression(expression)), Some(WorkflowStepInputSource.String(source))) => {
        val either =
          for {
            self <- lookupValue(source)
            pc = ParameterContext(inputValues, self)
            result <- expression.fold(EvaluateExpression).apply(pc, expressionLib).toEither
          } yield result

        either.toValidated
      }


      case (Some(StringOrExpression.Expression(expression)), Some(WorkflowStepInputSource.StringArray(array))) =>
        val either =
          for {
            self <- array.toList.traverse[ErrorOr, WomValue](lookupValue(_).toValidated).toEither
            pc = ParameterContext(inputValues, WomArray(self))
            result <- expression.fold(EvaluateExpression).apply(pc, expressionLib).toEither
          } yield result

        either.toValidated

      case (Some(StringOrExpression.Expression(expression)), None) =>
        expression.fold(EvaluateExpression).apply(ParameterContext(inputValues), expressionLib)

      case _ => s"Could not do evaluateValue(${input.valueFrom}, ${input.source}), most likely it has not been implemented yet".invalidNel
    }
  }

  //this is input, so not producing any output files
  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType) =
    Set.empty[WomFile].validNel

  override def inputs = graphInputs ++ input.source.toSet.flatMap{ inputSource: InputSource => inputSource match {
    case WorkflowStepInputSource.String(s) => Set(FullyQualifiedName(s).id)
    case WorkflowStepInputSource.StringArray(sa) => sa.map(FullyQualifiedName(_).id).toSet
  }}
}
