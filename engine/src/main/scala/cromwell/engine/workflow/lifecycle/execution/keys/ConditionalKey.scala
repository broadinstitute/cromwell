package cromwell.engine.workflow.lifecycle.execution.keys

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.ExecutionIndex.ExecutionIndex
import cromwell.core.logging.WorkflowLogger
import cromwell.core.{ExecutionStatus, JobKey}
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore.ValueKey
import cromwell.engine.workflow.lifecycle.execution.{WorkflowExecutionActorData, WorkflowExecutionDiff}
import wom.graph._
import wom.graph.expression.ExpressionNodeLike
import wom.values.{WomBoolean, WomValue}

/**
  * Represents a conditional node in the execution store.
  * Runnable when the associated expression (represented by an expression node in the graph) is done.
  */
private [execution] case class ConditionalKey(node: ConditionalNode, index: ExecutionIndex) extends JobKey {
  override val tag = node.localName
  override val attempt = 1

  lazy val collectors = node.conditionalOutputPorts map {
    ConditionalCollectorKey(_, index)
  }

  /**
    * Creates ExecutionStore entries for each of the scoped children.
    *
    * @return ExecutionStore of scattered children.
    */
  def populate(bypassed: Boolean): Map[JobKey, ExecutionStatus.Value] = {
    val conditionalKeys = node.innerGraph.nodes.flatMap({ node => keyify(node) })

    val finalStatus = if (bypassed) ExecutionStatus.NotStarted else ExecutionStatus.Bypassed
    (conditionalKeys ++ collectors).map({ _ -> finalStatus }).toMap
  }

  /**
    * Make a JobKey for all of the contained scopes.
    */
  private def keyify(node: GraphNode): Option[JobKey] = node match {
    case call: CommandCallNode => Option(BackendJobDescriptorKey(call, index, 1))
    case call: WorkflowCallNode => Option(SubWorkflowKey(call, index, 1))
    case expression: ExpressionNodeLike => Option(ExpressionKey(expression, index))
    case conditional: ConditionalNode => Option(ConditionalKey(conditional, index))
    case scatter: ScatterNode if index.isEmpty => Option(ScatterKey(scatter))
    case _: GraphInputNode => None
    case _: PortBasedGraphOutputNode => None
    case _: ScatterNode =>
      throw new UnsupportedOperationException("Nested Scatters are not supported (yet) ... but you might try a sub workflow to achieve the same effect!")
    case e =>
      throw new UnsupportedOperationException(s"Scope ${e.getClass.getName} is not supported in an If block.")
  }

  def processRunnable(data: WorkflowExecutionActorData, workflowLogger: WorkflowLogger): ErrorOr[WorkflowExecutionDiff] = {
    // This is the output port from the conditional's 'condition' input:
    val conditionOutputPort = node.conditionExpression.singleOutputPort
    data.valueStore.get(conditionOutputPort, index) match {
      case Some(b: WomBoolean) =>
        val conditionalStatus = if (b.value) ExecutionStatus.Done else ExecutionStatus.Bypassed

        val valueStoreAdditions: Map[ValueKey, WomValue] = if (!b.value) {
          node.conditionalOutputPorts.map(op => ValueKey(op, index) -> op.womType.none).toMap
        } else Map.empty

        WorkflowExecutionDiff(
          executionStoreChanges = populate(b.value) + (this -> conditionalStatus),
          valueStoreAdditions = valueStoreAdditions).validNel
      case Some(v: WomValue) =>
        s"'if' condition ${node.conditionExpression.womExpression.sourceString} must evaluate to a boolean but instead got ${v.womType.stableName}".invalidNel
      case None =>
        s"Could not find the boolean value for conditional $tag. Missing boolean should have come from expression ${node.conditionExpression.womExpression.sourceString}".invalidNel
    }
  }
}
