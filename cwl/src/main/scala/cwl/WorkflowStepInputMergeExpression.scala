package cwl

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
                                                  // cats doesn't have NonEmptyMap (yet https://github.com/typelevel/cats/pull/2141/)
                                                  // This is an ugly way to guarantee this class is only instantiated with at least one mapping
                                                  stepInputMappingHead: (String, OutputPort),
                                                  stepInputMappings: Map[String, OutputPort],
                                                  override val expressionLib: ExpressionLib) extends CwlWomExpression {
  private val allStepInputMappings = stepInputMappings + stepInputMappingHead

  override def sourceString: String = s"${input.id}-Merge-Expression"
  override def inputs: Set[String] = allStepInputMappings.keySet
  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
    def lookupValue(key: String): ErrorOr[WomValue] =
      inputValues.
        get(key).
        toValidNel(s"source value $key not found in input values ${inputValues.mkString("\n")}.  Graph Inputs were ${allStepInputMappings.keySet.mkString("\n")}")

    def validateSources(sources: List[String]): ErrorOr[List[WomValue]] =
      sources.
        traverse[ErrorOr, WomValue] {
        lookupValue
      }

    (input.source.map(_.fold(StringOrStringArrayToStringList)), input.effectiveLinkMerge) match {

      //When we have a single source, simply look it up
      case (Some(List(source)), LinkMergeMethod.MergeNested) =>
        lookupValue(source)

      //When we have several sources, validate they are all present and provide them as a nested array
      case (Some(sources), LinkMergeMethod.MergeNested) =>


        validateSources(sources).map(WomArray.apply)

      case (Some(sources), LinkMergeMethod.MergeFlattened) =>

        val validatedSourceValues: Checked[List[WomValue]] =
          validateSources(sources).toEither

        def flatten: WomValue => List[WomValue] = {
          case WomArray(_, value) => value.toList
          case WomOptionalValue(_, Some(value)) => flatten(value)
          case other => List(other)
        }

        //This is the meat of "merge_flattened," where we find arrays and concatenate them to form one array
        val flattenedValidatedSourceValues: Checked[List[WomValue]] = validatedSourceValues.map(_.flatMap {
          flatten
        })

        flattenedValidatedSourceValues.map(list => WomArray(list)).toValidated

      case (Some(List(id)), _) =>
        lookupValue(id)
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
