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
import wom.graph._
import wom.graph.expression._
import wom.types._

object WomBundleToFileElement {
  def convert(a: WomBundle): FileElement = {

    val tasks: Iterable[CallableTaskDefinition] = a.allCallables.values.filterByType[CallableTaskDefinition]
    val workflows: Iterable[WorkflowDefinition] = a.allCallables.values.filterByType[WorkflowDefinition]

    FileElement(
      Seq(),
      Seq(),
      workflows.map(WorkflowDefinitionToWorkflowDefinitionElement.convert).toSeq,
      tasks.map(CallableTaskDefinitionToTaskDefinitionElement.convert).toSeq
    )
  }
}

object CallableTaskDefinitionToTaskDefinitionElement {
  def convert(a: CallableTaskDefinition): TaskDefinitionElement = {
    TaskDefinitionElement(
      a.name,
      None,
      Seq(),
      None,
      CommandSectionElement(Seq()),
      None,
      MapToMetaSectionElement.convert(a.meta),
      MapToParameterMetaSectionElement.convert(a.parameterMeta)
    )
  }
}

object WorkflowDefinitionToWorkflowDefinitionElement {
  def convert(a: WorkflowDefinition): WorkflowDefinitionElement = {
    WorkflowDefinitionElement(
      a.name,
      None,
      a.graph.nodes.map(GraphNodeToWorkflowGraphElement.convert),
      None,
      MapToMetaSectionElement.convert(a.meta),
      MapToParameterMetaSectionElement.convert(a.parameterMeta)
    )
  }
}

// TODO: eliminate nested matches - they are impossible to read
object GraphNodeToWorkflowGraphElement {
  def convert(a: GraphNode): WorkflowGraphElement = {
    a match {
      case a: CallNode =>
        CallNodeToCallElement.convert(a)
      case a: ConditionalNode =>
        IfElement(
          conditionExpression = ExpressionNodeToExpressionElement.convert(a.conditionExpression),
          graphElements = Seq() //a.innerGraph.nodes
        )
      case a: ExpressionNodeLike =>
        a match {
          case a: ExpressionNode =>
            IntermediateValueDeclarationElement(
              typeElement = WomTypeToTypeElement.convert(a.womType),
              name = a.identifier.localName.value,
              expression = ExpressionNodeToExpressionElement.convert(a))
          case _: ExpressionCallNode => ???
        }
      case a: GraphNodeWithSingleOutputPort =>
        GraphNodeWithSingleOutputPortToWorkflowGraphElement.convert(a)
      case a: GraphOutputNode =>
        GraphOutputNodeToWorkflowGraphElement.convert(a)
      // a.scatterCollectionExpressionNodes.head.womExpression.sourceString
      // a.scatterCollectionExpressionNodes.head.womExpression.asInstanceOf[WdlWomExpression].wdlExpression
      case a: ScatterNode =>
        ScatterElement(
          scatterName = a.identifier.localName.value,
          scatterExpression = StringLiteral("wasd"),
          scatterVariableName = a.inputPorts.toList.head.name,
          graphElements = a.innerGraph.nodes.toList.map(GraphNodeToWorkflowGraphElement.convert)
        )
    }
  }
}
//object ExpressionNodeLikeToCallElement {
//  def convert(a: ExpressionNodeLike): CallElement = {
//    a match {
//      ExpressionCallNode
//    }
//  }
//}

object GraphOutputNodeToWorkflowGraphElement {
  def convert(a: GraphOutputNode): WorkflowGraphElement = {
    a match {
      case a: PortBasedGraphOutputNode =>
        OutputDeclarationElement(
          WomTypeToTypeElement.convert(a.womType),
          a.identifier.localName.value,
          StringLiteral("bogus") // TODO
        )
      case a: ExpressionBasedGraphOutputNode =>
        OutputDeclarationElement(
          WomTypeToTypeElement.convert(a.womType),
          a.identifier.localName.value,
          WomExpressionToExpressionElement.convert(a.womExpression)
        )
    }
  }
}

object GraphNodeWithSingleOutputPortToWorkflowGraphElement {
  def convert(a: GraphNodeWithSingleOutputPort): WorkflowGraphElement = {
    a match {
      case a: GraphInputNode =>
        InputDeclarationElement(
          WomTypeToTypeElement.convert(a.womType),
          a.identifier.localName.value,
          None
        )
      case a: ExpressionNode =>
        IntermediateValueDeclarationElement(
          WomTypeToTypeElement.convert(a.womType),
          a.identifier.localName.value,
          WomExpressionToExpressionElement.convert(a.womExpression)
        )
    }
  }
}

object WomTypeToTypeElement {
  def convert(a: WomType): TypeElement = {
    a match {
      case a: WomArrayType =>
        if (a.guaranteedNonEmpty)
          NonEmptyTypeElement(ArrayTypeElement(WomTypeToTypeElement.convert(a.memberType)))
        else
          ArrayTypeElement(WomTypeToTypeElement.convert(a.memberType))
      case _: WomCoproductType => ???
      case _: WomFileType => PrimitiveTypeElement(WomSingleFileType) // TODO: questionable assumption
      case a: WomMapType => MapTypeElement(WomTypeToTypeElement.convert(a.keyType), WomTypeToTypeElement.convert(a.valueType))
      case _: WomNothingType.type => ???
      case _: WomObjectType.type => ObjectTypeElement
      case a: WomOptionalType => OptionalTypeElement(WomTypeToTypeElement.convert(a.memberType))
      case a: WomPairType => PairTypeElement(WomTypeToTypeElement.convert(a.leftType), WomTypeToTypeElement.convert(a.rightType))
      case a: WomPrimitiveType => WomPrimitiveTypeToPrimitiveTypeElement.convert(a)
    }
  }
}

object WomPrimitiveTypeToPrimitiveTypeElement {
  def convert(a: WomPrimitiveType): PrimitiveTypeElement = PrimitiveTypeElement(a)
}

// TODO: Possibly the wrong conversion
object ExpressionNodeToExpressionElement {
  def convert(a: ExpressionNode): ExpressionElement = {
    WomExpressionToExpressionElement.convert(a.womExpression)
  }
}

object WomExpressionToExpressionElement {
  def convert(a: WomExpression): ExpressionElement = {
    a match {
      case _: WdlomWomExpression => ???
      //      case _: WdlWomExpression TODO: cannot import / does not make sense to have? Yet shows up.
      case _: ValueAsAnExpression => ???
      case _: InputLookupExpression => ???
      case _: PlainAnonymousExpressionNode => ???
      case _: TaskCallInputExpressionNode => ???
      case _: ExposedExpressionNode => ???
      //      case _ => throw new Exception("Unknown type")
      case _ => StringLiteral("todo")
    }
  }
}

object CallNodeToCallElement {
  def convert(a: CallNode): CallElement = {
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

// <<< Below this line, believed to be fully working >>>

object MapToMetaSectionElement {
  def convert(a: Map[String, String]): Option[MetaSectionElement] = {
    if (a.nonEmpty)
      Some(MetaSectionElement(a map { case (key, value) =>
        key -> MetaValueElementString(value) // draft-2: strings only
      }))
    else
      None
  }
}

object MapToParameterMetaSectionElement {
  def convert(a: Map[String, String]): Option[ParameterMetaSectionElement] = {
    if (a.nonEmpty)
      Some(ParameterMetaSectionElement(a map { case (key, value) =>
        key -> MetaValueElementString(value) // draft-2: strings only
      }))
    else
      None
  }
}
