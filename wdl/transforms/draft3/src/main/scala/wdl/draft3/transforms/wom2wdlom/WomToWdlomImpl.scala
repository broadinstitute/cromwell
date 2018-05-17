package wdl.draft3.transforms.wom2wdlom

import wdl.model.draft3.elements._
import wom.executable.WomBundle
import common.collections.EnhancedCollections.EnhancedTraversableLike
import wdl.draft2.model.WdlWomExpression
import wdl.draft3.transforms.wdlom2wom.expression.WdlomWomExpression
import wdl.model.draft3.elements.ExpressionElement.{ExpressionLiteralElement, StringLiteral}
import wdl.model.draft3.elements.MetaValueElement.MetaValueElementString
import wom.RuntimeAttributes
import wom.callable.{Callable, CallableTaskDefinition, WorkflowDefinition}
import wom.expression.{InputLookupExpression, ValueAsAnExpression, WomExpression}
import wom.graph.CallNode.InputDefinitionPointer
import wom.graph._
import wom.graph.expression._
import wom.types._

case object UnrepresentableException extends Exception("Value has no representation in the destination format (WDL)")

object WomToWdlomImpl {

  import WomToWdlom.ops._

  implicit val graphOutputNodeToWorkflowGraphElement: WomToWdlom[GraphOutputNode, WorkflowGraphElement] = {
    case a: PortBasedGraphOutputNode =>
      OutputDeclarationElement(
        a.womType.convert,
        a.identifier.localName.value,
        StringLiteral("bogus") // TODO
      )
    case a: ExpressionBasedGraphOutputNode =>
      OutputDeclarationElement(
        a.womType.convert,
        a.identifier.localName.value,
        a.womExpression.convert
      )
  }

  implicit val womBundleToFileElement: WomToWdlom[WomBundle, FileElement] = (a: WomBundle) => {
    val tasks: Iterable[CallableTaskDefinition] = a.allCallables.values.filterByType[CallableTaskDefinition]
    val workflows: Iterable[WorkflowDefinition] = a.allCallables.values.filterByType[WorkflowDefinition]

    FileElement(
      Seq(),
      Seq(),
      workflows.map(_.convert).toSeq,
      tasks.map(_.convert).toSeq
    )
  }

  implicit val mapToMetaSectionElement: WomToWdlom[Map[String, String], Option[MetaSectionElement]] = (a: Map[String, String]) => {
    if (a.nonEmpty)
      Some(MetaSectionElement(a map { case (key, value) =>
        key -> MetaValueElementString(value) // draft-2: strings only
      }))
    else
      None
  }

  implicit val mapToParameterMetaSectionElement: WomToWdlom[Map[String, String], Option[ParameterMetaSectionElement]] = (a: Map[String, String]) => {
    if (a.nonEmpty)
      Some(ParameterMetaSectionElement(a map { case (key, value) =>
        key -> MetaValueElementString(value) // draft-2: strings only
      }))
    else
      None
  }

  implicit val runtimeAttributesToRuntimeAttributesSectionElement: WomToWdlom[RuntimeAttributes, Option[RuntimeAttributesSectionElement]] = (a: RuntimeAttributes) => {
    def tupleToKvPair(tuple: (String, WomExpression)): ExpressionElement.KvPair =
      ExpressionElement.KvPair(tuple._1, tuple._2.convert)

    val kvPairs = (a.attributes map tupleToKvPair).toVector

    if (kvPairs.nonEmpty)
      Some(RuntimeAttributesSectionElement(kvPairs))
    else
      None
  }

  implicit val callableTaskDefinitionToTaskDefinitionElement: WomToWdlom[CallableTaskDefinition, TaskDefinitionElement] = (a: CallableTaskDefinition) => {
    TaskDefinitionElement(
      a.name,
      Some(InputsSectionElement(Seq())),
      Seq(),
      Some(OutputsSectionElement(Seq())),
      CommandSectionElement(Seq()),
      a.runtimeAttributes.convert,
      mapToMetaSectionElement.convert(a.meta), // TODO: why do these require explicit notation?
      mapToParameterMetaSectionElement.convert(a.parameterMeta)
    )
  }

  implicit val workflowDefinitionToWorkflowDefinitionElement: WomToWdlom[WorkflowDefinition, WorkflowDefinitionElement] = (a: WorkflowDefinition) => {
    WorkflowDefinitionElement(
      a.name,
      Some(InputsSectionElement(Seq())),
      a.graph.nodes.map(_.convert),
      Some(OutputsSectionElement(Seq())),
      mapToMetaSectionElement.convert(a.meta),
      mapToParameterMetaSectionElement.convert(a.parameterMeta)
    )
  }

  implicit val expressionNodeLikeToWorkflowGraphElement: WomToWdlom[ExpressionNodeLike, WorkflowGraphElement] = {
    case a: ExpressionNode =>
      IntermediateValueDeclarationElement(
        typeElement = a.womType.convert,
        name = a.identifier.localName.value,
        expression = a.convert)
    case _: ExpressionCallNode => ???
  }

  implicit val graphNodeToWorkflowGraphElement: WomToWdlom[GraphNode, WorkflowGraphElement] = {
    case a: CallNode =>
      a.convert
    case a: ConditionalNode =>
      IfElement(
        conditionExpression = a.conditionExpression.convert,
        graphElements = Seq() //a.innerGraph.nodes
      )
    case a: ExpressionNodeLike =>
      a.convert
    case a: GraphNodeWithSingleOutputPort =>
      a.convert
    case a: GraphOutputNode =>
      a.convert
    // a.scatterCollectionExpressionNodes.head.womExpression.sourceString
    // a.scatterCollectionExpressionNodes.head.womExpression.asInstanceOf[WdlWomExpression].wdlExpression
    case a: ScatterNode =>
      ScatterElement(
        scatterName = a.identifier.localName.value,
        scatterExpression = StringLiteral("wasd"),
        scatterVariableName = a.inputPorts.toList.head.name,
        graphElements = a.innerGraph.nodes.toList.map(graphNodeToWorkflowGraphElement.convert)
      )
  }

  implicit val graphNodeWithSingleOutputPortToWorkflowGraphElement: WomToWdlom[GraphNodeWithSingleOutputPort, WorkflowGraphElement] = {
    case a: GraphInputNode =>
      InputDeclarationElement(
        a.womType.convert,
        a.identifier.localName.value,
        None
      )
    case a: ExpressionNode =>
      IntermediateValueDeclarationElement(
        a.womType.convert,
        a.identifier.localName.value,
        a.womExpression.convert
      )
  }

  implicit val womTypeToTypeElement: WomToWdlom[WomType, TypeElement] = {
    case a: WomArrayType =>
      if (a.guaranteedNonEmpty)
        NonEmptyTypeElement(ArrayTypeElement(womTypeToTypeElement.convert(a.memberType)))
      else
        ArrayTypeElement(womTypeToTypeElement.convert(a.memberType))
    case _: WomCoproductType => throw UnrepresentableException
    case _: WomFileType => PrimitiveTypeElement(WomSingleFileType) // TODO: questionable assumption
    case a: WomMapType => MapTypeElement(womTypeToTypeElement.convert(a.keyType), womTypeToTypeElement.convert(a.valueType))
    case _: WomNothingType.type => throw UnrepresentableException
    case _: WomObjectType.type => ObjectTypeElement
    case a: WomOptionalType => OptionalTypeElement(womTypeToTypeElement.convert(a.memberType))
    case a: WomPairType => PairTypeElement(womTypeToTypeElement.convert(a.leftType), womTypeToTypeElement.convert(a.rightType))
    case a: WomPrimitiveType => a.convert
  }

  implicit val womPrimitiveTypeToPrimitiveTypeElement: WomToWdlom[WomPrimitiveType, PrimitiveTypeElement] = (a: WomPrimitiveType) => PrimitiveTypeElement(a)

  // TODO: Possibly the wrong conversion
  implicit val expressionNodeToExpressionElement: WomToWdlom[ExpressionNode, ExpressionElement] = (a: ExpressionNode) => {
    a.womExpression.convert
  }

  implicit val womExpressionToExpressionElement: WomToWdlom[WomExpression, ExpressionElement] = (a: WomExpression) => {
    a match {
      case _: WdlomWomExpression => ???
      case a: WdlWomExpression => ExpressionLiteralElement(a.sourceString)
      case _: ValueAsAnExpression => ???
      case _: InputLookupExpression => ???
      case _: PlainAnonymousExpressionNode => ???
      case _: TaskCallInputExpressionNode => ???
      case _: ExposedExpressionNode => ???
    }
  }

  implicit val callNodeToCallElement: WomToWdlom[CallNode, CallElement] = (a: CallNode) => {
    a.inputDefinitionMappings map { case (defn: Callable.InputDefinition, defnPtr: InputDefinitionPointer) =>
      (defn.localName, defnPtr)
    }
    a match {
      case a: ExpressionCallNode =>
        // We need to make list of inputs as kv pairs
        // An InputDefinition has a valueMapper, which seems like a callable function to obtain the input value?
        //        a.callable.inputs map { input =>
        //          input.valueMapper(NoIoFunctionSet)(???)
        //        }
        CallElement(
          a.callable.name,
          None,
          None
        )
      case _: CommandCallNode =>
        CallElement(
          a.callable.name,
          None,
          None
        )
      case _: WorkflowCallNode =>
        CallElement(
          a.callable.name,
          None,
          None
        )
    }
  }
}
