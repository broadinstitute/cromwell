package wdl.draft3.transforms.wom2wdlom

import wdl.model.draft3.elements.{WorkflowGraphElement, _}
import wom.executable.WomBundle
import common.collections.EnhancedCollections.EnhancedTraversableLike
import wdl.draft3.transforms.wdlom2wom.expression.WdlomWomExpression
import wdl.model.draft3.elements.ExpressionElement.StringLiteral
import wdl.model.draft3.elements.MetaValueElement.MetaValueElementString
import wom.callable.{Callable, CallableTaskDefinition, WorkflowDefinition}
import wom.expression.{InputLookupExpression, ValueAsAnExpression, WomExpression}
import wom.graph.CallNode.InputDefinitionPointer
import wom.graph.{GraphNodeWithSingleOutputPort, _}
import wom.graph.expression._
import wom.types._

case object UnrepresentableException extends Exception("Value has no representation in the destination format (WDL)")

object ConvertImpl {

  import Convert.ops._

  implicit val graphOutputNodeToWorkflowGraphElement: Convert[GraphOutputNode, WorkflowGraphElement] = {
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

  implicit val womBundleToFileElement: Convert[WomBundle, FileElement] = (a: WomBundle) => {
    val tasks: Iterable[CallableTaskDefinition] = a.allCallables.values.filterByType[CallableTaskDefinition]
    val workflows: Iterable[WorkflowDefinition] = a.allCallables.values.filterByType[WorkflowDefinition]

    FileElement(
      Seq(),
      Seq(),
      workflows.map(_.convert).toSeq,
      tasks.map(_.convert).toSeq
    )
  }

  implicit val mapToMetaSectionElement: Convert[Map[String, String], Option[MetaSectionElement]] = (a: Map[String, String]) => {
    if (a.nonEmpty)
      Some(MetaSectionElement(a map { case (key, value) =>
        key -> MetaValueElementString(value) // draft-2: strings only
      }))
    else
      None
  }

  implicit val mapToParameterMetaSectionElement: Convert[Map[String, String], Option[ParameterMetaSectionElement]] = (a: Map[String, String]) => {
    if (a.nonEmpty)
      Some(ParameterMetaSectionElement(a map { case (key, value) =>
        key -> MetaValueElementString(value) // draft-2: strings only
      }))
    else
      None
  }

  implicit val callableTaskDefinitionToTaskDefinitionElement: Convert[CallableTaskDefinition, TaskDefinitionElement] = (a: CallableTaskDefinition) => {
    TaskDefinitionElement(
      a.name,
      None,
      Seq(),
      None,
      CommandSectionElement(Seq()),
      None,
      mapToMetaSectionElement.convert(a.meta), // TODO: why do these require explicit notation?
      mapToParameterMetaSectionElement.convert(a.parameterMeta)
    )
  }

  implicit val workflowDefinitionToWorkflowDefinitionElement: Convert[WorkflowDefinition, WorkflowDefinitionElement] = (a: WorkflowDefinition) => {
    WorkflowDefinitionElement(
      a.name,
      None,
      a.graph.nodes.map(_.convert),
      None,
      mapToMetaSectionElement.convert(a.meta),
      mapToParameterMetaSectionElement.convert(a.parameterMeta)
    )
  }

  implicit val expressionNodeLikeToWorkflowGraphElement: Convert[ExpressionNodeLike, WorkflowGraphElement] = {
    case a: ExpressionNode =>
      IntermediateValueDeclarationElement(
        typeElement = a.womType.convert,
        name = a.identifier.localName.value,
        expression = a.convert)
    case _: ExpressionCallNode => ???
  }

  implicit val graphNodeToWorkflowGraphElement: Convert[GraphNode, WorkflowGraphElement] = {
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

  implicit val graphNodeWithSingleOutputPortToWorkflowGraphElement: Convert[GraphNodeWithSingleOutputPort, WorkflowGraphElement] = {
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

  implicit val womTypeToTypeElement: Convert[WomType, TypeElement] = {
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

  implicit val womPrimitiveTypeToPrimitiveTypeElement: Convert[WomPrimitiveType, PrimitiveTypeElement] = (a: WomPrimitiveType) => PrimitiveTypeElement(a)

  // TODO: Possibly the wrong conversion
  implicit val expressionNodeToExpressionElement: Convert[ExpressionNode, ExpressionElement] = (a: ExpressionNode) => {
    a.womExpression.convert
  }

  implicit val womExpressionToExpressionElement: Convert[WomExpression, ExpressionElement] = {
    case _: WdlomWomExpression => ???
    //      case a: wdl.draft2.model.WdlWomExpression
    //      case _: WdlWomExpression TODO: cannot import / does not make sense to have? Yet shows up at runtime.
    case _: ValueAsAnExpression => ???
    case _: InputLookupExpression => ???
    case _: PlainAnonymousExpressionNode => ???
    case _: TaskCallInputExpressionNode => ???
    case _: ExposedExpressionNode => ???
    //      case _ => throw new Exception("Unknown type")
    case _ => StringLiteral("todo")
  }

  implicit val callNodeToCallElement: Convert[CallNode, CallElement] = (a: CallNode) => {
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
