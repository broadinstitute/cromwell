package wdl.draft3.transforms.wom2wdlom

import wdl.model.draft3.elements._
import wom.executable.WomBundle
import common.collections.EnhancedCollections.EnhancedTraversableLike
import wdl.model.draft3.elements.ExpressionElement.StringLiteral
import wdl.model.draft3.elements.MetaValueElement.MetaValueElementString
import wom.callable.{CallableTaskDefinition, WorkflowDefinition}
import wom.graph._
import wom.graph.expression.ExpressionNodeLike

object WomBundleToFileElement {
  def convert(a: WomBundle): FileElement = {

    val tasks: Iterable[CallableTaskDefinition] = a.allCallables.values.filterByType[CallableTaskDefinition]
    val workflows: Iterable[WorkflowDefinition] = a.allCallables.values.filterByType[WorkflowDefinition]

    FileElement(
      Seq(),
      Seq(),
      workflows.map(WorkflowDefinitionToTaskDefinitionElement.convert).toSeq,
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

object WorkflowDefinitionToTaskDefinitionElement {
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

object GraphNodeToWorkflowGraphElement {
  def convert(a: GraphNode): WorkflowGraphElement = {
    a match {
      case a: CallNode =>
        CallNodeToCallElement.convert(a)
      case a: ConditionalNode =>
        ScatterElement(
          scatterName = a.identifier.localName.value,
          scatterExpression = StringLiteral("wasd"),
          scatterVariableName = "a",
          graphElements = Seq()
        )
      case a: ExpressionNodeLike =>
        ScatterElement(
          scatterName = a.identifier.localName.value,
          scatterExpression = StringLiteral("wasd"),
          scatterVariableName = "a",
          graphElements = Seq()
        )
      case a: GraphNodeWithSingleOutputPort =>
        ScatterElement(
          scatterName = a.identifier.localName.value,
          scatterExpression = StringLiteral("wasd"),
          scatterVariableName = "a",
          graphElements = Seq()
        )
      case a: GraphOutputNode =>
        ScatterElement(
          scatterName = a.identifier.localName.value,
          scatterExpression = StringLiteral("wasd"),
          scatterVariableName = "a",
          graphElements = Seq()
        )
      case a: ScatterNode =>
        ScatterElement(
          scatterName = a.identifier.localName.value,
          scatterExpression = StringLiteral("wasd"),
          scatterVariableName = "a",
          graphElements = Seq()
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

object CallNodeToCallElement {
  def convert(a: CallNode): CallElement = {
    a match {
      case a: ExpressionCallNode =>
        a.callable.inputs map { input =>
          input.valueMapper

        }
        CallElement(
          a.callable.name,
          None,
          None
        )
      case _: CommandCallNode => ???
      case _: WorkflowCallNode => ???
    }
  }
}

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
