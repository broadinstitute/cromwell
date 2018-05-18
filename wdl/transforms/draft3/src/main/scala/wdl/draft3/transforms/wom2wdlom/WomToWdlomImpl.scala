package wdl.draft3.transforms.wom2wdlom

import wdl.model.draft3.elements._
import wom.executable.WomBundle
import common.collections.EnhancedCollections.EnhancedTraversableLike
import wdl.draft2.model.WdlWomExpression
import wdl.draft3.transforms.wdlom2wom.expression.WdlomWomExpression
import wom.callable.Callable._
import wdl.model.draft3.elements.ExpressionElement.{ExpressionLiteralElement, StringLiteral}
import wdl.model.draft3.elements.MetaValueElement.MetaValueElementString
import wom.RuntimeAttributes
import wom.callable.Callable.OutputDefinition
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
        a.womType.toWdlom,
        a.identifier.localName.value,
        StringLiteral("bogus") // TODO
      )
    case a: ExpressionBasedGraphOutputNode =>
      OutputDeclarationElement(
        a.womType.toWdlom,
        a.identifier.localName.value,
        a.womExpression.toWdlom
      )
  }

  implicit val womBundleToFileElement: WomToWdlom[WomBundle, FileElement] = (a: WomBundle) => {
    val tasks: Iterable[CallableTaskDefinition] = a.allCallables.values.filterByType[CallableTaskDefinition]
    val workflows: Iterable[WorkflowDefinition] = a.allCallables.values.filterByType[WorkflowDefinition]

    FileElement(
      Seq(),
      Seq(),
      workflows.map(_.toWdlom).toSeq,
      tasks.map(_.toWdlom).toSeq
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
      ExpressionElement.KvPair(tuple._1, tuple._2.toWdlom)

    val kvPairs = (a.attributes map tupleToKvPair).toVector

    if (kvPairs.nonEmpty)
      Some(RuntimeAttributesSectionElement(kvPairs))
    else
      None
  }

  implicit val outputDefinitionToOutputDeclarationElement: WomToWdlom[OutputDefinition, OutputDeclarationElement] = (a: OutputDefinition) => {
    OutputDeclarationElement(a.womType.toWdlom, a.localName.value, a.expression.toWdlom)
  }

  implicit val inputDefinitionToInputDeclarationElement: WomToWdlom[InputDefinition, InputDeclarationElement] = {
    case a: RequiredInputDefinition => InputDeclarationElement(a.womType.toWdlom, a.localName.value, None)
    case a: InputDefinitionWithDefault => InputDeclarationElement(a.womType.toWdlom, a.localName.value, Some(a.default.toWdlom))
    case a: FixedInputDefinition => InputDeclarationElement(a.womType.toWdlom, a.localName.value, Some(a.default.toWdlom))
    case a: OptionalInputDefinition => InputDeclarationElement(a.womType.toWdlom, a.localName.value, None)
  }

  implicit val callableTaskDefinitionToTaskDefinitionElement: WomToWdlom[CallableTaskDefinition, TaskDefinitionElement] = (a: CallableTaskDefinition) => {
    // TODO: why does _.toWdlom not work here?
    val inputs = a.inputs.map(inputDefinitionToInputDeclarationElement.toWdlom)
    val outputs = a.outputs.map(outputDefinitionToOutputDeclarationElement.toWdlom)

    TaskDefinitionElement(
      a.name,
      if (inputs.nonEmpty) Some(InputsSectionElement(inputs)) else None,
      Seq(),
      if (outputs.nonEmpty) Some(OutputsSectionElement(outputs)) else None,
      CommandSectionElement(Seq()),
      a.runtimeAttributes.toWdlom,
      mapToMetaSectionElement.toWdlom(a.meta), // TODO: why do these require explicit notation?
      mapToParameterMetaSectionElement.toWdlom(a.parameterMeta)
    )
  }

  implicit val workflowDefinitionToWorkflowDefinitionElement: WomToWdlom[WorkflowDefinition, WorkflowDefinitionElement] = (a: WorkflowDefinition) => {
    val inputs = a.inputs.map(inputDefinitionToInputDeclarationElement.toWdlom)
    val outputs = a.outputs.map(outputDefinitionToOutputDeclarationElement.toWdlom)

    WorkflowDefinitionElement(
      a.name,
      if (inputs.nonEmpty) Some(InputsSectionElement(inputs)) else None,
      a.graph.nodes.map(_.toWdlom),
      if (outputs.nonEmpty) Some(OutputsSectionElement(outputs)) else None,
      mapToMetaSectionElement.toWdlom(a.meta),
      mapToParameterMetaSectionElement.toWdlom(a.parameterMeta)
    )
  }

  implicit val expressionNodeLikeToWorkflowGraphElement: WomToWdlom[ExpressionNodeLike, WorkflowGraphElement] = {
    case a: ExpressionNode =>
      IntermediateValueDeclarationElement(
        typeElement = a.womType.toWdlom,
        name = a.identifier.localName.value,
        expression = a.toWdlom)
    case _: ExpressionCallNode => ???
  }

  implicit val graphNodeToWorkflowGraphElement: WomToWdlom[GraphNode, WorkflowGraphElement] = {
    case a: CallNode =>
      a.toWdlom
    case a: ConditionalNode =>
      IfElement(
        conditionExpression = a.conditionExpression.toWdlom,
        graphElements = Seq() //a.innerGraph.nodes
      )
    case a: ExpressionNodeLike =>
      a.toWdlom
    case a: GraphNodeWithSingleOutputPort =>
      a.toWdlom
    case a: GraphOutputNode =>
      a.toWdlom
    // a.scatterCollectionExpressionNodes.head.womExpression.sourceString
    // a.scatterCollectionExpressionNodes.head.womExpression.asInstanceOf[WdlWomExpression].wdlExpression
    case a: ScatterNode =>
      ScatterElement(
        scatterName = a.identifier.localName.value,
        scatterExpression = StringLiteral("wasd"),
        scatterVariableName = a.inputPorts.toList.head.name,
        graphElements = a.innerGraph.nodes.toList.map(graphNodeToWorkflowGraphElement.toWdlom)
      )
  }

  implicit val graphNodeWithSingleOutputPortToWorkflowGraphElement: WomToWdlom[GraphNodeWithSingleOutputPort, WorkflowGraphElement] = {
    case a: GraphInputNode =>
      InputDeclarationElement(
        a.womType.toWdlom,
        a.identifier.localName.value,
        None
      )
    case a: ExpressionNode =>
      IntermediateValueDeclarationElement(
        a.womType.toWdlom,
        a.identifier.localName.value,
        a.womExpression.toWdlom
      )
  }

  implicit val womTypeToTypeElement: WomToWdlom[WomType, TypeElement] = {
    case a: WomArrayType =>
      if (a.guaranteedNonEmpty)
        NonEmptyTypeElement(ArrayTypeElement(womTypeToTypeElement.toWdlom(a.memberType)))
      else
        ArrayTypeElement(womTypeToTypeElement.toWdlom(a.memberType))
    case _: WomCoproductType => throw UnrepresentableException
    case _: WomFileType => PrimitiveTypeElement(WomSingleFileType) // TODO: questionable assumption
    case a: WomMapType => MapTypeElement(womTypeToTypeElement.toWdlom(a.keyType), womTypeToTypeElement.toWdlom(a.valueType))
    case _: WomNothingType.type => throw UnrepresentableException
    case _: WomObjectType.type => ObjectTypeElement
    case a: WomOptionalType => a.toWdlom
    case a: WomPairType => PairTypeElement(womTypeToTypeElement.toWdlom(a.leftType), womTypeToTypeElement.toWdlom(a.rightType))
    case a: WomPrimitiveType => a.toWdlom
  }

  implicit val womOptionalTypeToOptionalTypeElement: WomToWdlom[WomOptionalType, OptionalTypeElement] = (a: WomOptionalType) =>
    OptionalTypeElement(womTypeToTypeElement.toWdlom(a.memberType))

  implicit val womPrimitiveTypeToPrimitiveTypeElement: WomToWdlom[WomPrimitiveType, PrimitiveTypeElement] = (a: WomPrimitiveType) => PrimitiveTypeElement(a)

  // TODO: Possibly the wrong conversion
  implicit val expressionNodeToExpressionElement: WomToWdlom[ExpressionNode, ExpressionElement] = (a: ExpressionNode) => {
    a.womExpression.toWdlom
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
