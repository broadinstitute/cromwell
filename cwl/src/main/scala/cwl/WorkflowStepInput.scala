package cwl

import common.validation.ErrorOr.ErrorOr
import cwl.LinkMergeMethod.LinkMergeMethod
import cwl.WorkflowStepInput.InputSource
import cwl.command.ParentName
import wom.graph.GraphNodePort.OutputPort
import wom.graph.WomIdentifier
import wom.graph.expression.ExposedExpressionNode
import wom.types.WomType
import cats.syntax.option._
import cats.syntax.traverse._
import cats.instances.list._
import eu.timepit.refined._
import cats.syntax.either._
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
    * @param expectedType the type which the step's run is expecting
    * @return
    */
  def toExpressionNode(sourceMappings: Map[String, OutputPort],
                       outputTypeMap: Map[String, WomType],
                       expressionLib: ExpressionLib,
                       expectedType: Option[MyriadInputType],
                       isScattered: Boolean
                      )(implicit parentName: ParentName): ErrorOr[ExposedExpressionNode] = {

    val sources = source.toList.flatMap(_.fold(StringOrStringArrayToStringList)).map(FullyQualifiedName(_).id)

    val inputs = sourceMappings.keySet

    val outputTypeMapWithIDs = outputTypeMap.map {
      case (key, value) => FullyQualifiedName(key).id -> value
    }

    def lookupId(id: String): ErrorOr[WomType] =
      outputTypeMapWithIDs.
        get(id).
        toValidNel(s"couldn't find $id as derived from $source in map\n${outputTypeMapWithIDs.mkString("\n")}")

    (for {
      //lookup each of our source Ids, failing if any of them are missing
      inputTypes <- sources.traverse[ErrorOr, WomType](lookupId).toEither
      //we may have several sources, we make sure to have a type common to all of them
      inputType: WomType = expectedType.map(_.fold(MyriadInputTypeToWomType)) getOrElse WomType.homogeneousTypeFromTypes(inputTypes)
      womExpression = WorkflowStepInputExpression(this, inputType, inputs, expressionLib)
      identifier = WomIdentifier(id)
      ret <- ExposedExpressionNode.fromInputMapping(identifier, womExpression, inputType, sourceMappings).toEither
    } yield ret).toValidated
  }
}

object WorkflowStepInput {
  type InputSource = String :+: Array[String] :+: CNil
}