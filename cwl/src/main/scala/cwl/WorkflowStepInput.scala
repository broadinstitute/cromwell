package cwl

import cats.data.NonEmptyList
import common.validation.ErrorOr.ErrorOr
import cwl.LinkMergeMethod.LinkMergeMethod
import cwl.WorkflowStepInput.InputSource
import cwl.command.ParentName
import wom.graph.GraphNodePort.OutputPort
import wom.graph.WomIdentifier
import wom.graph.expression.ExposedExpressionNode
import wom.types.{WomArrayType, WomType}
import cats.syntax.traverse._
import cats.syntax.option._
import cats.instances.list._
import eu.timepit.refined._
import cats.syntax.either._
import common.Checked
import shapeless.{:+:, CNil, Witness}
import shapeless.syntax.singleton._
import cwl.LinkMergeMethod.LinkMergeMethod
import cwl.WorkflowStepInput.InputSource
import common.validation.ErrorOr.ErrorOr
import cwl.CommandLineTool.{CommandBindingSortingKey, SortKeyAndCommandPart}
import cwl.command.ParentName
import wom.types.WomType
import wom.graph.WomIdentifier
import wom.graph.GraphNodePort.OutputPort
import wom.graph.expression.ExposedExpressionNode
import wom.values.WomValue

case class WorkflowStepInput(
  id: String,
  source: Option[InputSource] = None,
  linkMerge: Option[LinkMergeMethod] = None,
  default: Option[CwlAny] = None,
  valueFrom: Option[StringOrExpression] = None) {

  /**
    *
    * @param sourceMappings The outputports to which this source refers
    * @param outputTypeMap The types of the output ports to which this source refers
    * @param matchingRunInputType This input matches an input declared in the workflowstep's "run".  This is that step's declared type
    * @return
    */
  def toExpressionNode(sourceMappings: Map[String, OutputPort],
                       outputTypeMap: Map[String, WomType],
                       expressionLib: ExpressionLib,
                       matchingRunInputType: Option[MyriadInputType],
                       isScattered: Boolean
                      )(implicit parentName: ParentName): ErrorOr[ExposedExpressionNode] = {

    val inputs = sourceMappings.keySet

    val identifier = WomIdentifier(id)

    val node =
      for {
        inputType <- WorkflowStepInput.determineType(this, outputTypeMap, matchingRunInputType, isScattered)
        womExpression = WorkflowStepInputExpression(this, inputType, inputs, expressionLib)
        node <- ExposedExpressionNode.fromInputMapping(identifier, womExpression, inputType, sourceMappings).toEither
      } yield node
    node.toValidated
  }
}

object WorkflowStepInput {
  type InputSource = String :+: Array[String] :+: CNil

  def determineType(stepInput: WorkflowStepInput,
                    outputTypeMap: Map[String, WomType],
                    expectedType: Option[MyriadInputType],
                    isScattered: Boolean)(implicit parentName: ParentName): Checked[WomType] = {

    val sources = stepInput.source.toList.flatMap(_.fold(StringOrStringArrayToStringList)).map(FullyQualifiedName(_).id)

    /*
    val outputTypeMapWithIDs = outputTypeMap.map {
      case (key, value) => FullyQualifiedName(key).id -> value
    }

    def lookupId(id: String): ErrorOr[WomType] =
      outputTypeMapWithIDs.
        get(id).
        toValidNel(s"couldn't find $id as derived from ${stepInput.source} in map\n${outputTypeMapWithIDs.mkString("\n")}")

    lazy val typeFromSources: Either[NonEmptyList[String], List[WomType]] = sources.traverse[ErrorOr, WomType](lookupId).toEither

    val typeFromRunInput: Option[WomType] = expectedType.map(_.fold(MyriadInputTypeToWomType))

    val sinkParameterIsArray: Boolean = expectedType.map(_ match {
      case MyriadInputType.InputArray(_) => true
      case _ => false
    }).getOrElse(false)
    */

    (sources.size,isScattered, expectedType) match {
      case (n, true, Some(tpe)) if (n > 0) => WomArrayType(tpe.fold(MyriadInputTypeToWomType)).asRight
    }
  }
}