package cromwell.engine.workflow.lifecycle.execution.keys

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.{ExecutionStatus, JobKey}
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore.ValueKey
import cromwell.engine.workflow.lifecycle.execution.{WorkflowExecutionActorData, WorkflowExecutionDiff}
import wom.graph._
import wom.graph.expression.{ExposedExpressionNode, ExpressionNode}
import wom.values.WomArray.WomArrayLike
import wom.values.WomValue

import scala.language.postfixOps

private [execution] case class ScatterKey(node: ScatterNode) extends JobKey {
  // When scatters are nested, this might become Some(_)
  override val index = None
  override val attempt = 1
  override val tag = node.localName

  def makeCollectors(count: Int): Set[ScatterCollectorKey] = (node.outputMapping.groupBy(_.outputToGather.source.graphNode) flatMap {
    case (_: CallNode | _: ExposedExpressionNode | _: ConditionalNode, scatterGatherPorts) => scatterGatherPorts.map(sgp => ScatterCollectorKey(sgp, count))
    case _ => Set.empty[ScatterCollectorKey]
  }).toSet

  /**
    * Creates a sub-ExecutionStore with Starting entries for each of the scoped children.
    *
    * @param count Number of ways to scatter the children.
    * @return ExecutionStore of scattered children.
    */
  def populate(count: Int): Map[JobKey, ExecutionStatus.Value] = {
    val shards = node.innerGraph.nodes flatMap { makeShards(_, count) }
    val collectors = makeCollectors(count)
    (shards ++ collectors) map { _ -> ExecutionStatus.NotStarted } toMap
  }

  private def makeShards(scope: GraphNode, count: Int): Seq[JobKey] = scope match {
    case call: TaskCallNode => (0 until count) map { i => BackendJobDescriptorKey(call, Option(i), 1) }
    case expression: ExpressionNode => (0 until count) map { i => ExpressionKey(expression, Option(i)) }
    case conditional: ConditionalNode => (0 until count) map { i => ConditionalKey(conditional, Option(i)) }
    case subworkflow: WorkflowCallNode => (0 until count) map { i => SubWorkflowKey(subworkflow, Option(i), 1) }
    case _: GraphInputNode => List.empty
    case _: PortBasedGraphOutputNode => List.empty
    case _: ScatterNode =>
      throw new UnsupportedOperationException("Nested Scatters are not supported (yet) ... but you might try a sub workflow to achieve the same effect!")
    case e =>
      throw new UnsupportedOperationException(s"Scope ${e.getClass.getName} is not supported.")
  }

  def processRunnable(data: WorkflowExecutionActorData): ErrorOr[WorkflowExecutionDiff] = {
    val collectionOutputPort = node.scatterCollectionExpressionNode.singleExpressionOutputPort

    data.valueStore.get(collectionOutputPort, None) map {
      case WomArrayLike(arrayLike) =>
        WorkflowExecutionDiff(
          // Add the new shards + collectors, and set the scatterVariable and scatterKey to Done
          executionStoreChanges = populate(arrayLike.value.size) ++ Map(
            ScatterVariableInputKey(node.scatterVariableInnerGraphInputNode, arrayLike) -> ExecutionStatus.Done,
            this -> ExecutionStatus.Done
          ),
          // Add scatter variable arrayLike to the value store
          valueStoreAdditions = Map(
            ValueKey(node.scatterVariableInnerGraphInputNode.singleOutputPort, None) -> arrayLike
          )
        ).validNel
      case v: WomValue =>
        s"Scatter collection ${node.scatterCollectionExpressionNode.womExpression.sourceString} must evaluate to an array but instead got ${v.womType.toDisplayString}".invalidNel
    } getOrElse {
      s"Could not find an array value for scatter $tag. Missing array should have come from expression ${node.scatterCollectionExpressionNode.womExpression.sourceString}".invalidNel
    }
  }
}
