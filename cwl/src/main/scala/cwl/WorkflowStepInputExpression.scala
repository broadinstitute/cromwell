package cwl

import cats.data.NonEmptyList
import cats.syntax.validated._
import cats.syntax.either._
import cats.syntax.traverse._
import cats.syntax.option._
import cats.instances.list._
import common.Checked
import common.validation.ErrorOr.ErrorOr
import cwl.WorkflowStepInput.InputSource
import cwl.command.ParentName
import mouse.all._
import wom.expression.IoFunctionSet
import wom.types._
import wom.values._

final case class WorkflowStepInputExpression(input: WorkflowStepInput,
                                             override val cwlExpressionType: WomType,
                                             graphInputs: Set[String],
                                             override val expressionLib: ExpressionLib)(implicit parentName: ParentName) extends CwlWomExpression {

  override def sourceString = input.toString

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet) = {

    def lookupValue(key: String): ErrorOr[WomValue] =
      inputValues.
        get(FullyQualifiedName(key).id).
        toValidNel(s"source value $key not found in input values ${inputValues.mkString("\n")}.  Graph Inputs were ${graphInputs.mkString("\n")}")

    (input.valueFrom, input.source.map(_.fold(StringOrStringArrayToStringList)), input.effectiveLinkMerge) match {
        //When we have a single source, simply look it up
      case (None, Some(List(source)), LinkMergeMethod.MergeNested) =>
        lookupValue(source)

      //When we have a single source, validate they are all present and provide them as a nested array of maps
      case (None, Some(sources), LinkMergeMethod.MergeNested) =>

        def convertTupleToWomMap: (String, WomValue) => WomMap = {
          case (key, value) => WomMap(WomMapType(WomStringType, value.womType), Map(WomString(key) -> value))
        }

        val validatedSourceValues: ErrorOr[List[WomMap]] =
          sources.
            traverse[ErrorOr, WomMap]{
              s =>
                lookupValue(s).
                  // bind each wom value to its original id so it can be referred to by name
                  map(FullyQualifiedName(s).id -> _).

                  // convert each tuple into a map consisting of a single element
                  map(convertTupleToWomMap.tupled)
            }

        validatedSourceValues.map(WomArray.apply)

      case (None, Some(sources), LinkMergeMethod.MergeFlattened) =>

        val validatedSourceValues: Checked[List[WomValue]] = sources.
          traverse[ErrorOr, WomValue](s => lookupValue(s)).
          toEither

        def flatten: WomValue => List[WomValue] = {
          womValue =>
            womValue match {
              case WomArray(_, value) => value.toList
              case WomOptionalValue(_, Some(value)) => flatten(value)
              case other => List(other)
            }
        }
        //This is the meat of "merge_flattened," where we find arrays and concatenate them to form one array
        val flattenedValidatedSourceValues: Checked[List[WomValue]] = validatedSourceValues.map(list => list.flatMap{ womValue =>
          flatten(womValue)
        })

        flattenedValidatedSourceValues.map(list => WomArray(list)).toValidated

      case (None, Some(List(id)), _) =>
        lookupValue(id)

      // If valueFrom is a constant string value, use this as the value for this input parameter.
      // TODO: need to handle case where this is a parameter reference, it currently looks like a String to us!
      case (Some(StringOrExpression.String(value)), None, _) => WomString(value).validNel

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
      case (Some(StringOrExpression.Expression(expression)), Some(sources), _) =>

        //used to determine the value of "self" as expected by CWL Spec
        def selfValue(in: List[WomValue]) = in match {
          case single :: Nil => single
          case multiple => WomArray(multiple)
        }

        val evaluatedExpression: Checked[WomValue] =
          for {
            sourceValueList <- sources.traverse[ErrorOr, WomValue](lookupValue).toEither
            //value is either a womArray or not depending on how many items are in the source list
            self = selfValue(sourceValueList)
            pc = ParameterContext(inputValues, self)
            result <- expression.fold(EvaluateExpression).apply(pc, expressionLib).toEither
          } yield result

        evaluatedExpression.toValidated

      case (Some(StringOrExpression.Expression(expression)), None, _) =>
        expression.fold(EvaluateExpression).apply(ParameterContext(inputValues), expressionLib)

      case other => s"Could not do evaluateValue $other, most likely it has not been implemented yet".invalidNel
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


