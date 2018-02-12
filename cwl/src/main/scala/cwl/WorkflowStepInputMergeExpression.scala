package cwl

import cats.data.NonEmptyList
import cats.instances.list._
import cats.syntax.either._
import cats.syntax.option._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.Checked
import common.validation.ErrorOr.ErrorOr
import wom.expression.IoFunctionSet
import wom.graph.GraphNodePort.OutputPort
import wom.types.WomType
import wom.values.{WomArray, WomFile, WomOptionalValue, WomValue}

final case class WorkflowStepInputMergeExpression(input: WorkflowStepInput,
                                                  cwlExpressionType: WomType,
                                                  stepInputMappingHead: NonEmptyList[(String, OutputPort)],
                                                  override val expressionLib: ExpressionLib) extends CwlWomExpression {

  private val allStepInputMappings = stepInputMappingHead.toList
  private val allStepInputSources = allStepInputMappings.map(_._1)

  override def sourceString: String = s"${input.id}-Merge-Expression"
  override def inputs: Set[String] = allStepInputSources.toSet

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
    def lookupValue(key: String): ErrorOr[WomValue] =
      inputValues.
        get(key).
        toValidNel(s"source value $key not found in input values ${inputValues.mkString("\n")}.  Graph Inputs were ${allStepInputSources.mkString("\n")}")

    def validateSources(sources: List[String]): ErrorOr[List[WomValue]] =
      sources.
        traverse[ErrorOr, WomValue](lookupValue)

    (allStepInputSources, input.effectiveLinkMerge) match {

      //When we have a single source, simply look it up
      case (List(source), LinkMergeMethod.MergeNested) => lookupValue(source)

      //When we have several sources, validate they are all present and provide them as a nested array
      case (sources, LinkMergeMethod.MergeNested) => validateSources(sources).map(WomArray.apply)

      case (sources, LinkMergeMethod.MergeFlattened) =>
        val validatedSourceValues: Checked[List[WomValue]] =
          validateSources(sources).toEither

        def flatten: WomValue => List[WomValue] = {
          case WomArray(_, value) => value.toList
          case WomOptionalValue(_, Some(value)) => flatten(value)
          case other => List(other)
        }

        //This is the meat of "merge_flattened," where we find arrays and concatenate them to form one array
        val flattenedValidatedSourceValues: Checked[List[WomValue]] = validatedSourceValues.map(_.flatMap(flatten))

        flattenedValidatedSourceValues.map(list => WomArray(list)).toValidated

      case (List(id), _) => lookupValue(id)
    }
  }

  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] = {
    if (allStepInputMappings.size > 1) {
      // TODO add MultipleInputFeatureRequirement logic in here
      "MultipleInputFeatureRequirement not supported yet".invalidNel
    } else {
      val (inputName, _) = allStepInputMappings.head
      inputTypes(inputName).collectAsSeq({
        case file: WomFile => file
      }).toSet.validNel
    }
  }
}
