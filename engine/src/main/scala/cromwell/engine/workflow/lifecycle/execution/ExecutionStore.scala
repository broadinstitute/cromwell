package cromwell.engine.workflow.lifecycle.execution

import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.ExecutionIndex.ExecutionIndex
import cromwell.core.ExecutionStatus._
import cromwell.core.JobKey
import cromwell.engine.workflow.lifecycle.execution.ExecutionStore.{RunnableScopes, _}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.{apply => _}
import cromwell.engine.workflow.lifecycle.execution.keys._
import lenthall.collections.Table
import wom.callable.WorkflowDefinition
import wom.graph.GraphNodePort.{ConditionalOutputPort, OutputPort, ScatterGathererPort}
import wom.graph._
import wom.graph.expression.ExpressionNode

object ExecutionStore {
  implicit class EnhancedOutputPort(val outputPort: OutputPort) extends AnyVal {
    def executionNode: GraphNode = outputPort match {
      case scatter: ScatterGathererPort => scatter.outputToGather.executionNode
      case conditional: ConditionalOutputPort => conditional.outputToExpose.executionNode
      case other => other.graphNode.executionNode
    }
  }

  implicit class EnhancedGraphNode(val graphNode: GraphNode) extends AnyVal {
    def executionNode: GraphNode = graphNode match {
      case pbgon: PortBasedGraphOutputNode => pbgon.source.graphNode
      case other => other
    }
  }
  
  case class RunnableScopes(scopes: List[JobKey], truncated: Boolean)

  def empty = ExecutionStore(Map.empty[JobKey, ExecutionStatus], hasNewRunnables = false)

  def apply(workflow: WorkflowDefinition): ExecutionStore = {
    // Keys that are added in a NotStarted Status
    val notStartedKeys = workflow.innerGraph.nodes collect {
      case call: TaskCallNode => List(BackendJobDescriptorKey(call, None, 1))
      case expression: ExpressionNode => List(ExpressionKey(expression, None))
      case scatterNode: ScatterNode => List(ScatterKey(scatterNode))
      case conditionalNode: ConditionalNode => List(ConditionalKey(conditionalNode, None))
    }

    // There are potentially resolved workflow inputs that are default WomExpressions.
    // For now assume that those are call inputs that will be evaluated in the CallPreparation.
    // If they are actually workflow declarations then we would need to add them to the ExecutionStore so they can be evaluated.
    // In that case we would want InstantiatedExpressions so we can create an InstantiatedExpressionNode and add a DeclarationKey
    new ExecutionStore(notStartedKeys.flatten.map(_ -> NotStarted).toMap, notStartedKeys.nonEmpty)
  }

  val MaxJobsToStartPerTick = 1000
}

final case class ExecutionStore(private val statusStore: Map[JobKey, ExecutionStatus], hasNewRunnables: Boolean) {
  // View of the statusStore more suited for lookup based on status
  lazy val store: Map[ExecutionStatus, List[JobKey]] = statusStore.groupBy(_._2).mapValues(_.keys.toList)

  /*
    * Create 2 Tables, one for keys in done status and one for keys in terminal status.
    * A Table is nothing more than a Map[R, Map[C, V]], see Table trait for more details
    * In this case, rows are GraphNodes, columns are ExecutionIndexes, and values are JobKeys
    * This allows for quick lookup of all shards for a node, as well as accessing a specific key with a 
    * (node, index) pair
   */
  lazy val (doneKeys, terminalKeys) = {
    def toTableEntry(key: JobKey) = (key.node, key.index, key)
    store.foldLeft((Table.empty[GraphNode, ExecutionIndex, JobKey], Table.empty[GraphNode, ExecutionIndex, JobKey]))({
      case ((done, terminal), (status, keys))  =>
        lazy val newMapEntries = keys map toTableEntry
        val newDone = if (status.isDoneOrBypassed) done.addAll(newMapEntries) else done
        val newTerminal = if (status.isTerminal) terminal.addAll(newMapEntries) else terminal

        newDone -> newTerminal
    })
  }

  private def keysWithStatus(status: ExecutionStatus) = store.getOrElse(status, List.empty)

  def isInBypassedConditional(jobKey: JobKey): Boolean = keysWithStatus(Bypassed).exists {
    case conditional: ConditionalKey if conditional.node.innerGraph.nodes.contains(jobKey.node)
      && conditional.index.equals(jobKey.index) => true
    case _ => false
  }

  def hasActiveJob: Boolean = {
    def upstreamFailed(node: GraphNode): Boolean = node.upstreamAncestry exists hasFailedNode

    keysWithStatus(QueuedInCromwell).nonEmpty ||
      keysWithStatus(Starting).nonEmpty ||
      keysWithStatus(Running).nonEmpty ||
      keysWithStatus(NotStarted).exists(jobKey => !upstreamFailed(jobKey.node))
  }

  def jobStatus(jobKey: JobKey): Option[ExecutionStatus] = statusStore.get(jobKey)

  def startedJobs: List[BackendJobDescriptorKey] = {
    store.filterNot({ case (s, _) => s == NotStarted}).values.toList.flatten collect {
      case k: BackendJobDescriptorKey => k
    }
  }

  private def hasFailedNode(node: GraphNode): Boolean = keysWithStatus(Failed).exists(_.node == node)

  def hasFailedJob: Boolean = keysWithStatus(Failed).nonEmpty

  override def toString: String = store.map { case (j, s) => s"$j -> $s" } mkString System.lineSeparator()

  def add(values: Map[JobKey, ExecutionStatus]): ExecutionStore = {
    this.copy(statusStore = statusStore ++ values, hasNewRunnables = hasNewRunnables || values.values.exists(_.isTerminalOrRetryable))
  }

  /**
    * Returns the list of jobs ready to be run, along with a Boolean indicating whether or not the list has been truncated.
    * The size of the list will be MaxJobsToStartPerTick at most. If more jobs where found runnable, the boolean will be true, otherwise false.
    */
  def runnableScopes: RunnableScopes = {
    val readyToStart = keysWithStatus(NotStarted).toStream filter arePrerequisitesDone
    // Compute the first ExecutionStore.MaxJobsToStartPerTick + 1 runnable scopes
    val scopesToStartPlusOne = readyToStart.take(ExecutionStore.MaxJobsToStartPerTick + 1).toList
    // Only take the first ExecutionStore.MaxJobsToStartPerTick from the above list.
    // Use the fact that we took one more to determine whether or not we truncated the result.
    RunnableScopes(scopesToStartPlusOne.take(ExecutionStore.MaxJobsToStartPerTick), scopesToStartPlusOne.size > ExecutionStore.MaxJobsToStartPerTick)
  }

  private def arePrerequisitesDone(key: JobKey): Boolean = {
    def calculateUpstreamDone(upstream: GraphNode, index: ExecutionIndex) = upstream match {
        // FIXME this won't hold conditionals in scatters for example
        // OuterGraphInputNode signals that an input comes from outside the graph.
        // Depending on whether or not this input is outside of a scatter graph will change the index which we need to look at
      case outerNode: OuterGraphInputNode => doneKeys.contains(outerNode, None)
      case _: CallNode | _: ScatterNode | _: ExpressionNode | _: ConditionalNode => doneKeys.contains(upstream, index)
      case _ => true
    }

    def upstreamAreDone = key.node.upstreamPorts forall {
      // The collector is at index None, so if this is a scatter gather port ignore the key index
      case upstreamPort: ScatterGathererPort => calculateUpstreamDone(upstreamPort.executionNode, None)
      case upstreamPort => calculateUpstreamDone(upstreamPort.executionNode, key.index)
    }

    val runnable = key match {
      case scatterCollector: ScatterCollectorKey => terminalKeys.row(scatterCollector.node).size == scatterCollector.scatterWidth
      case conditionalCollector: ConditionalCollectorKey => terminalKeys.contains(conditionalCollector.node.executionNode, key.index)
      case _ => upstreamAreDone
    }

    runnable
  }
}
