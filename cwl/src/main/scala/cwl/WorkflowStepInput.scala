package cwl

import cats.data.NonEmptyList
import cats.syntax.either._
import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import cwl.InputParameter.DefaultToWomValuePoly
import cwl.LinkMergeMethod.LinkMergeMethod
import cwl.WorkflowStepInput.InputSource
import cwl.command.ParentName
import shapeless.{:+:, CNil}
import wom.expression.ValueAsAnExpression
import wom.graph.GraphNodePort.OutputPort
import wom.graph.WomIdentifier
import wom.graph.expression.{AnonymousExpressionNode, ExposedExpressionNode, ExpressionNode, PlainAnonymousExpressionNode}
import wom.types._

case class WorkflowStepInput(
                              id: String,
                              source: Option[InputSource] = None,
                              linkMerge: Option[LinkMergeMethod] = None,
                              default: Option[CwlAny] = None,
                              valueFrom: Option[StringOrExpression] = None) {

  def parsedId(implicit parentName: ParentName) = FullyQualifiedName(id).id

  def toExpressionNode(valueFromExpression: StringOrExpression,
                       runInputExpectedType: Option[cwl.MyriadInputType],
                       isScattered: Boolean,
                       sourceMappings:Map[String, OutputPort],
                       outputTypeMap: Map[String, WomType],
                       expressionLib: ExpressionLib,
                       schemaDefRequirement: SchemaDefRequirement
                      )(implicit parentName: ParentName): ErrorOr[ExpressionNode] = {
    val inputs = sourceMappings.keySet
    val upstreamMergeType = outputTypeMap.get(parsedId)

    (for {
      // we may have several sources, we make sure to have a type common to all of them.
      // In the case where there's no input source, we currently wrap the valueFrom value in a WomString (see WorkflowStepInputExpression)
      inputType <- WorkflowStepInput.determineValueFromType(upstreamMergeType, runInputExpectedType, isScattered, schemaDefRequirement)
      womExpression = WorkflowStepInputExpression(parsedId, valueFromExpression, inputType, inputs, expressionLib)
      identifier = WomIdentifier(id).combine("expression")
      ret <- ExposedExpressionNode.fromInputMapping(identifier, womExpression, inputType, sourceMappings).toEither
    } yield ret).toValidated
  }

  /**
    *
    * @param sourceMappings The outputports to which this source refers
    * @param matchingRunInputType This input matches an input declared in the workflowstep's "run".  This is that step's declared type
    * @return
    */
  def toMergeNode(sourceMappings: Map[String, OutputPort],
                  expressionLib: ExpressionLib,
                  matchingRunInputType: Option[MyriadInputType],
                  isScattered: Boolean,
                  schemaDefRequirement: SchemaDefRequirement
                 ): Option[ErrorOr[ExpressionNode]] = {

    val identifier = WomIdentifier(id).combine("merge")
    val mapType = sourceMappings.map({ case (k, v) => k -> v.womType })

    val maybeMatchingRunInputWomType: Option[WomType] = matchingRunInputType.map(_.fold(MyriadInputTypeToWomType).apply(schemaDefRequirement))

    def makeNode(head: (String, OutputPort), tail: List[(String, OutputPort)]) = for {
      inputType <- determineMergeType(mapType, maybeMatchingRunInputWomType)
      womExpression = WorkflowStepInputMergeExpression(this, inputType, NonEmptyList.of(head, tail: _*), expressionLib)
      node <- AnonymousExpressionNode.fromInputMapping(identifier, womExpression, sourceMappings, PlainAnonymousExpressionNode.apply).toEither
    } yield node

    val matchingRunInputWomType: WomType = maybeMatchingRunInputWomType.getOrElse(WomAnyType)
    lazy val defaultValue = default.map { _.fold(DefaultToWomValuePoly).apply(matchingRunInputWomType) }

    sourceMappings.toList match {
      case head :: tail => Option(makeNode(head, tail).toValidated)
      case Nil =>
        defaultValue map { _ map { d =>
            PlainAnonymousExpressionNode(identifier, ValueAsAnExpression(d), matchingRunInputWomType, Map.empty)
        } }
    }
  }

  def determineMergeType(sources: Map[String, WomType], expectedTypeAsWom: Option[WomType]): Checked[WomType] = {
    WorkflowStepInput.determineMergeType(sources, linkMerge, expectedTypeAsWom, default.isDefined)
  }

  lazy val sources: List[String] = source.toList.flatMap(_.fold(StringOrStringArrayToStringList))

  lazy val effectiveLinkMerge: LinkMergeMethod = linkMerge.getOrElse(LinkMergeMethod.MergeNested)
}

object WorkflowStepInput {
  type InputSource = String :+: Array[String] :+: CNil

  implicit class EnhancedStepInputMap[A](val map: Map[WorkflowStepInput, A]) extends AnyVal {
    def asIdentifierMap(implicit parentName: ParentName): Map[String, A] = {
      map.map({ case (stepInput, value) => stepInput.parsedId -> value })
    }
  }

  def determineValueFromType(mergedSourcesType: Option[WomType],
                             expectedType: Option[MyriadInputType],
                             isScattered: Boolean,
                             schemaDefRequirement: SchemaDefRequirement): Checked[WomType] = {
    val expectedTypeAsWom: Option[WomType] = expectedType.map(_.fold(MyriadInputTypeToWomType).apply(schemaDefRequirement))

    expectedTypeAsWom.getOrElse(WomStringType).asRight
  }

  def determineMergeType(sources: Map[String, WomType],
                         linkMerge: Option[LinkMergeMethod],
                         expectedTypeAsWom: Option[WomType],
                         hasDefault: Boolean): Checked[WomType] = {

    (sources.toList, expectedTypeAsWom, linkMerge) match {
      // If there is a single source and no explicit merge method, use the type of the source (unpacked if it's optional
      // and there's a default since the default would be used if there's no supplied value).
      case (List((_, WomOptionalType(sourceType))), _, None) if hasDefault => sourceType.asRight
      case (List((_, opt @ WomOptionalType(_))), _, None) => opt.asRight
      case (List((_, sourceType)), _, None) => sourceType.asRight
      // If there is a single source and merge nested is specified, wrap it into an array
      case (List((_, sourceType)), _, Some(LinkMergeMethod.MergeNested)) => WomArrayType(sourceType).asRight
      // If there are multiple sources and either no merge method or merge nested, find the closest common type.
      // TODO: This is technically not enough, cwltool supports sources with totally different types, creating an array with multiple types
      // Maybe WomCoproduct can help
      case (_, _, Some(LinkMergeMethod.MergeNested) | None) => WomArrayType(WomType.homogeneousTypeFromTypes(sources.values)).asRight

      //If sink parameter is an array and merge_flattened is used, must validate input & output types are equivalent before proceeding
      case (_, Some(arrayType @ WomArrayType(itemType)), Some(LinkMergeMethod.MergeFlattened)) if typesToItemMatch(sources.values, itemType) => arrayType.asRight
      // If sink parameter is not an array and merge flattened is used, validate that the sources types matche the sink type
      case (_, Some(targetType), Some(LinkMergeMethod.MergeFlattened)) if typesToItemMatch(sources.values, targetType) => WomArrayType(targetType).asRight
      // If the types are not compatible, fail
      case (_, Some(targetType), Some(LinkMergeMethod.MergeFlattened)) =>
        s"could not verify that types $sources and the items type of the run's InputArraySchema $targetType were compatible".invalidNelCheck

      //We don't have type information from the run input so we gather up the sources and try to determine a common type amongst them.
      case _ => WomType.homogeneousTypeFromTypes(sources.values).asRight
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
