package cwl

import cats.data.NonEmptyList
import cats.instances.list._
import cats.syntax.either._
import cats.syntax.option._
import cats.syntax.traverse._
import common.Checked
import common.validation.ErrorOr.ErrorOr
import cwl.LinkMergeMethod.LinkMergeMethod
import cwl.WorkflowStepInput.InputSource
import cwl.command.ParentName
import mouse.all._
import shapeless.{:+:, CNil}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.WomIdentifier
import wom.graph.expression.{AnonymousExpressionNode, ExposedExpressionNode, ExpressionNode, PlainAnonymousExpressionNode}
import wom.types.{WomArrayType, WomStringType, WomType}
import wom.values.WomValue

case class WorkflowStepInput(
  id: String,
  source: Option[InputSource] = None,
  linkMerge: Option[LinkMergeMethod] = None,
  default: Option[CwlAny] = None,
  valueFrom: Option[StringOrExpression] = None) {

  def parsedId(implicit parentName: ParentName) = FullyQualifiedName(id).id

  def toExpressionNode(valueFromExpression: StringOrExpression,
                       sourceMappings:Map[String, OutputPort],
                       outputTypeMap: Map[String, WomType],
                       expressionLib: ExpressionLib
                      )(implicit parentName: ParentName): ErrorOr[ExpressionNode] = {
    val inputs = sourceMappings.keySet

    def lookupId(id: String): ErrorOr[WomType] =
      outputTypeMap.
        get(id).
        toValidNel(s"couldn't find $id as derived from $source in map\n${outputTypeMap.mkString("\n")}")
    

    (for {
      //lookup each of our source Ids, failing if any of them are missing
      inputTypes <- inputs.toList.traverse[ErrorOr, WomType](lookupId).toEither
      // we may have several sources, we make sure to have a type common to all of them.
      // In the case where there's no input source, we currently wrap the valueFrom value in a WomString (see WorkflowStepInputExpression)
      inputType = if (inputTypes.isEmpty) WomStringType else WomType.homogeneousTypeFromTypes(inputTypes)
      womExpression = WorkflowStepInputExpression(parsedId, valueFromExpression, inputType, inputs, expressionLib)
      identifier = WomIdentifier(id).combine("expression")
      ret <- ExposedExpressionNode.fromInputMapping(identifier, womExpression, inputType, sourceMappings).toEither
    } yield ret).toValidated
  }

  /**
    *
    * @param sourceMappings The outputports to which this source refers
    * @param outputTypeMap The types of the output ports to which this source refers
    * @param matchingRunInputType This input matches an input declared in the workflowstep's "run".  This is that step's declared type
    * @return
    */
  def toMergeNode(sourceMappings: Map[String, OutputPort],
                       outputTypeMap: Map[String, WomType],
                       expressionLib: ExpressionLib,
                       matchingRunInputType: Option[MyriadInputType],
                       isScattered: Boolean
                      )(implicit parentName: ParentName): Option[ErrorOr[ExpressionNode]] = {
    
    val identifier = WomIdentifier(id).combine("merge")

    def makeNode(head: (String, OutputPort), tail: List[(String, OutputPort)]) = for {
        inputType <- WorkflowStepInput.determineType(this, outputTypeMap, matchingRunInputType, isScattered)
        womExpression = WorkflowStepInputMergeExpression(this, inputType, head, tail.toMap, expressionLib)
        node <- AnonymousExpressionNode.fromInputMapping(identifier, womExpression, sourceMappings, PlainAnonymousExpressionNode.apply).toEither
      } yield node

    sourceMappings.toList match {
      case Nil => None
      case head :: tail => Option(makeNode(head, tail).toValidated)
    }
  }

  lazy val sources: List[String] = source.toList.flatMap(_.fold(StringOrStringArrayToStringList))

  lazy val effectiveLinkMerge: LinkMergeMethod = linkMerge.getOrElse(LinkMergeMethod.MergeNested)

  def sourceValues(inputValues: Map[String, WomValue])(implicit pn: ParentName): Checked[Map[String, WomValue]]  = {

    def lookupValue(key: String): Checked[WomValue] =
      inputValues.
        get(FullyQualifiedName(key).id).
        toRight(s"source value $key not found in input values ${inputValues.mkString("\n")}." |> NonEmptyList.one)


    sources.
      traverse[ErrorOr, (String, WomValue)](s => lookupValue(s).toValidated.map(FullyQualifiedName(s).id -> _)).
      toEither.
      map(_.toMap)
  }

  def validatedSourceTypes(typeMap: WomTypeMap)(implicit pn: ParentName): Checked[WomTypeMap] = {
    def lookupValue(key: String): Checked[WomType] =
      typeMap.
        get(FullyQualifiedName(key).id).
        toRight(s"source value $key not found in type map ${typeMap.mkString("\n")}." |> NonEmptyList.one)

    sources.
      traverse[ErrorOr, (String, WomType)](s => lookupValue(s).toValidated.map(FullyQualifiedName(s).id -> _)).
      toEither.
      map(_.toMap)
  }

}

object WorkflowStepInput {
  type InputSource = String :+: Array[String] :+: CNil

  implicit class EnhancedStepInputMap[A](val map: Map[WorkflowStepInput, A]) extends AnyVal {
    def asIdentifierMap(implicit parentName: ParentName): Map[String, A] = {
      map.map({ case (stepInput, value) => stepInput.parsedId -> value })
    }
  }

  def determineType(stepInput: WorkflowStepInput,
                    outputTypeMap: Map[String, WomType],
                    expectedType: Option[MyriadInputType],
                    isScattered: Boolean)(implicit parentName: ParentName): Checked[WomType] = {

    val expectedTypeAsWom: Option[WomType] = expectedType.map(_.fold(MyriadInputTypeToWomType))
    (isScattered, expectedTypeAsWom, stepInput.effectiveLinkMerge) match {

      //If scattering over this variable, we expect an array of the sink type
      case (true, Some(tpe), _) => WomArrayType(tpe).asRight

      //If sink parameter is an array, we must frame the input as an array
      case (false, Some(array: WomArrayType), LinkMergeMethod.MergeNested) =>
        array.asRight

      //If sink parameter is an array and merge_flattened is used, must validate input & output types are equivalent before proceeding
      case (false, Some(array@WomArrayType(tpe)), LinkMergeMethod.MergeFlattened)  =>
        //Collect the upstream outputs and their types
        stepInput.validatedSourceTypes(outputTypeMap) flatMap {
          case map  if typesToItemMatch(map.values, tpe) =>  array.asRight
          case map => (s"could not verify that types $map and the items type of the run's InputArraySchema $tpe were compatible" |> NonEmptyList.one).asLeft
        }

      //We don't need to alter the type in the base case (i.e. not a scatter variable and no array type in the run input)
      case (false, Some(tpe), _) => tpe.asRight

      //We don't have type information from the run input so we gather up the sources and try to determine a common type amongst them.
      case _ => stepInput.validatedSourceTypes(outputTypeMap).map(_.values).map(WomType.homogeneousTypeFromTypes)
    }
  }

  def typesToItemMatch(lst: Iterable[WomType], target: WomType): Boolean = {
    val effectiveInputType = WomType.homogeneousTypeFromTypes(lst)

    typeToItemMatch(effectiveInputType, target)
  }

  def typeToItemMatch(upstream: WomType, downstream: WomType): Boolean = {
    upstream match {
      case WomType.RecursiveType(innerType) => typeToItemMatch(innerType, downstream)
      case other => other == downstream
    }
  }
}
