package cromwell.engine.workflow.lifecycle.execution.keys

import cats.syntax.either._
import cats.syntax.validated._
import common.Checked
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.{ExecutionStatus, JobKey}
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore.ValueKey
import cromwell.engine.workflow.lifecycle.execution.{WorkflowExecutionActorData, WorkflowExecutionDiff}
import wom.graph.ScatterNode.{ScatterCollectionFunction, ScatterVariableAndValue}
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

  /**
    * Creates a sub-ExecutionStore with Starting entries for each of the scoped children.
    *
    * @param count Number of ways to scatter the children.
    * @return ExecutionStore of scattered children.
    */
  def populate(count: Int, scatterCollectionFunction: ScatterCollectionFunction): Map[JobKey, ExecutionStatus.Value] = {
    val shards = node.innerGraph.nodes flatMap { makeShards(_, count) }

    /*
      * There's quite a few levels of indirection, but this is roughly what it looks like:
      * ScatterNode
      * |    * CallNode A
      * |      |- OutputPort A
      * |             ^
      * |             |
      * |    * OutputToGather A (PortBasedGraphOutputNode) -> one per output per shard
      * |             ^
      * |             |
      * |- ScatterGathererPort -> CollectedValue (array) of all shards for a given output port
      * 
      * This represents a ScatterNode containing one call (CallNodeA) which itself exposes one output port (OutputPort A).
      * There could be multiple output ports for the same CallNode and multiple CallNodes (or Conditional or Declaration Nodes etc...) inside the scatter.
      * The scatter node contains a PortBasedGraphOutputNode for each output port it wants to expose to the outside.
      * For each of those PortBasedGraphOutputNode, there is an associated ScatterGathererPort (they're the outputMapping of the scatter node).
      * 
      * We could make a collector for each scatter gather port (as they all represent an output of the scatter).
      * However in order to limit the number of those collectors (which are treated like normal nodes)
      * we can optimize a little by making only one collector if the node at the other end of the scatter gather port
      * is a TaskCallNode. This is because we expect all the output ports of a TaskCallNode to be available at the same time,
      * indeed we get all the outputs for a task call at once (when the call is complete). If this assumption is no longer
      * true and we get partial output values while a task is running this optimization would need to be removed.
      * We make that same assumption for WorkflowCallNodes (sub workflows) but this could be changed as well as it's not
      * unlikely that we'd want workflow outputs to be available as soon as possible.
      * 
      * 
      * In order to do this optimization, we group the scatter gather ports by node they point to
      * (not the PortBasedGraphOutputNode that is just a proxy, but the source node). If that node is a CallNode
      * we create only one collector for all those scatter gather ports. 
     */
    val collectors: Set[ScatterCollectorKey] = (node.outputMapping.groupBy(_.outputToGather.source.graphNode) flatMap {
      case (_: CallNode | _: ExposedExpressionNode | _: ConditionalNode, scatterGatherPorts) => scatterGatherPorts.map(sgp => ScatterCollectorKey(sgp, count, scatterCollectionFunction))
      case _ => Set.empty[ScatterCollectorKey]
    }).toSet

    (shards ++ collectors) map { _ -> ExecutionStatus.NotStarted } toMap
  }

  private def makeShards(scope: GraphNode, count: Int): Seq[JobKey] = scope match {
    case call: TaskCallNode => (0 until count) map { i => BackendJobDescriptorKey(call, Option(i), 1) }
    case expression: ExpressionNode => (0 until count) map { i => ExpressionKey(expression, Option(i)) }
    case conditional: ConditionalNode => (0 until count) map { i => ConditionalKey(conditional, Option(i)) }
    case subworkflow: WorkflowCallNode => (0 until count) map { i => SubWorkflowKey(subworkflow, Option(i), 1) }
    case _: GraphInputNode => List.empty
    case _: PortBasedGraphOutputNode => List.empty
    //        case call: WdlWorkflowCall => (0 until count) map { i => SubWorkflowKey(call, Option(i), 1) }
    case _: ScatterNode =>
      throw new UnsupportedOperationException("Nested Scatters are not supported (yet) ... but you might try a sub workflow to achieve the same effect!")
    case e =>
      throw new UnsupportedOperationException(s"Scope ${e.getClass.getName} is not supported.")
  }

  def processRunnable(data: WorkflowExecutionActorData): ErrorOr[WorkflowExecutionDiff] = {
    import cats.instances.list._
    import cats.syntax.traverse._

    def getScatterArray(scatterVariableNode: ScatterVariableNode): ErrorOr[ScatterVariableAndValue] = {
      val expressionNode = scatterVariableNode.scatterExpressionNode
      data.valueStore.get(expressionNode.singleExpressionOutputPort, None) map {
        case WomArrayLike(arrayLike) => ScatterVariableAndValue(scatterVariableNode, arrayLike).validNel
        case v: WomValue =>
          s"Scatter collection ${expressionNode.womExpression.sourceString} must evaluate to an array but instead got ${v.womType.toDisplayString}".invalidNel
      } getOrElse {
        s"Could not find an array value for scatter $tag. Missing array should have come from expression ${expressionNode.womExpression.sourceString}".invalidNel
      }
    }

    // The scatter iteration nodes mapped to their value retrieved from the value store
    val scatterArraysValuesCheck: Checked[List[ScatterVariableAndValue]] = node
      // Get all the iteration nodes (there will be as many as variables we're scattering over)
      .scatterVariableNodes
      // Retrieve the values of the collection nodes value from the ValueStore
      .traverse[ErrorOr, ScatterVariableAndValue](getScatterArray)
      // Convert to either so we can flatMap later
      .toEither

    // Execution changes (for execution store and value store) generated by the scatter iteration nodes
    def buildExecutionChanges(scatterVariableAndValues: List[ScatterVariableAndValue]) = {
      val (executionStoreChanges, valueStoreChanges) = scatterVariableAndValues.map({
        case ScatterVariableAndValue(scatterVariableNode, arrayValue) =>
          val executionStoreChange = ScatterVariableInputKey(scatterVariableNode, arrayValue) -> ExecutionStatus.Done
          val valueStoreChange = ValueKey(scatterVariableNode.singleOutputPort, None) -> arrayValue

          executionStoreChange -> valueStoreChange
      }).unzip

      executionStoreChanges.toMap -> valueStoreChanges.toMap
    }

    (for {
      arrays <- scatterArraysValuesCheck
      scatterSize <- node.scatterProcessingFunction(arrays)
      (executionStoreChanges, valueStoreChanges) = buildExecutionChanges(arrays)
    } yield WorkflowExecutionDiff(
      executionStoreChanges = populate(scatterSize, node.scatterCollectionFunctionBuilder(arrays.map(_.arrayValue.size))) ++ executionStoreChanges ++ Map(this -> ExecutionStatus.Done),
      valueStoreAdditions = valueStoreChanges
    )).toValidated
  }
}
