package cromwell.engine.workflow.lifecycle.execution

import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.ExecutionStatus._
import cromwell.core.{CallKey, JobKey}
import cromwell.engine.workflow.lifecycle.execution.ExecutionStore.FqnIndex
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.{apply => _, _}
import wdl4s._

import scala.language.postfixOps

object ExecutionStore {
  type FqnIndex = (String, Option[Int])

  def empty = ExecutionStore(Map.empty[JobKey, ExecutionStatus], hasNewRunnables = false)

  def apply(workflow: Workflow, workflowCoercedInputs: WorkflowCoercedInputs) = {
    // Only add direct children to the store, the rest is dynamically created when necessary
    val keys = workflow.children map {
      case call: TaskCall => Option(BackendJobDescriptorKey(call, None, 1))
      case call: WorkflowCall => Option(SubWorkflowKey(call, None, 1))
      case scatter: Scatter => Option(ScatterKey(scatter))
      case conditional: If => Option(ConditionalKey(conditional, None))
      case declaration: Declaration => Option(DeclarationKey(declaration, None, workflowCoercedInputs))
      case _ => None
    }

    new ExecutionStore(keys.flatten.map(_ -> NotStarted).toMap, keys.nonEmpty)
  }
  
  val MaxJobsToStartPerTick = 1000
}

final case class ExecutionStore(private val statusStore: Map[JobKey, ExecutionStatus], hasNewRunnables: Boolean) {

  // View of the statusStore more suited for lookup based on status
  lazy val store: Map[ExecutionStatus, List[JobKey]] = statusStore.groupBy(_._2).mapValues(_.keys.toList)
  // Takes only keys that are done, and creates a map such that they're indexed by fqn and index
  // This allows for quicker lookup (by hash) instead of traversing the whole list and yields
  // significant improvements at large scale (run ExecutionStoreBenchmark)
  lazy val doneKeys: Map[FqnIndex, JobKey] = store.filterKeys(_.isDoneOrBypassed).values.flatten.map { key =>
    (key.scope.fullyQualifiedName, key.index) -> key
  } toMap

  private def keysWithStatus(status: ExecutionStatus) = store.getOrElse(status, List.empty)

  def isBypassedConditional(jobKey: JobKey, conditional: If): Boolean = {
    keysWithStatus(Bypassed).exists {
      case key: ConditionalKey =>
        key.scope.fullyQualifiedName.equals(conditional.fullyQualifiedName) &&
          key.index.equals(jobKey.index)
      case _ => false
    }
  }

  def hasActiveJob: Boolean = {
    def upstreamFailed(scope: Scope): Boolean = scope match {
      case node: GraphNode => node.upstreamAncestry exists hasFailedScope
    }

    keysWithStatus(QueuedInCromwell).nonEmpty ||
      keysWithStatus(Starting).nonEmpty ||
      keysWithStatus(Running).nonEmpty ||
      keysWithStatus(NotStarted).exists(jobKey => !upstreamFailed(jobKey.scope))
  }

  def jobStatus(jobKey: JobKey): Option[ExecutionStatus] = statusStore.get(jobKey)

  def startedJobs: List[BackendJobDescriptorKey] = {
    store.filterNot({ case (s, _) => s == NotStarted}).values.toList.flatten collect {
      case k: BackendJobDescriptorKey => k
    }
  }

  private def hasFailedScope(s: GraphNode): Boolean = keysWithStatus(Failed).exists(_.scope == s)

  def hasFailedJob: Boolean = keysWithStatus(Failed).nonEmpty

  override def toString = store.map { case (j, s) => s"$j -> $s" } mkString System.lineSeparator()

  def add(values: Map[JobKey, ExecutionStatus]) = {
    this.copy(statusStore = statusStore ++ values, hasNewRunnables = hasNewRunnables || values.values.exists(_.isTerminal))
  }

  /**
    * Returns the list of jobs ready to be run, along with a Boolean indicating whether or not the list has been truncated.
    * The size of the list will be MaxJobsToStartPerTick at most. If more jobs where found runnable, the boolean will be true, otherwise false.
    */
  def runnableScopes: (List[JobKey], Boolean) = {
    val readyToStart = keysWithStatus(NotStarted).toStream filter arePrerequisitesDone
    (readyToStart.take(ExecutionStore.MaxJobsToStartPerTick).toList, readyToStart.size > ExecutionStore.MaxJobsToStartPerTick)
  }

  def findCompletedShardsForOutput(key: CollectorKey): List[JobKey] = doneKeys.values.toList collect {
    case k @ (_: CallKey | _:DynamicDeclarationKey) if k.scope == key.scope && k.isShard => k
  }

  private def emulateShardEntries(key: CollectorKey): Set[FqnIndex] = {
    (0 until key.scatterWidth).toSet map { i: Int => key.scope match {
      case c: Call => c.fullyQualifiedName -> Option(i)
      case d: Declaration => d.fullyQualifiedName -> Option(i)
      case _ => throw new RuntimeException("Don't collect that.")
    }}
  }

  private def arePrerequisitesDone(key: JobKey): Boolean = {
    lazy val upstreamAreDone = key.scope.upstream forall {
      case n @ (_: Call | _: Scatter | _: Declaration) => upstreamIsDone(key, n)
      case _ => true
    }

    val shardEntriesForCollectorAreDone: Boolean = key match {
      case collector: CollectorKey => emulateShardEntries(collector).diff(doneKeys.keys.toSet).isEmpty
      case _ => true
    }

    shardEntriesForCollectorAreDone && upstreamAreDone
  }

  private def upstreamIsDone(entry: JobKey, prerequisiteScope: Scope): Boolean = {
    prerequisiteScope.closestCommonAncestor(entry.scope) match {
      /*
        * If this entry refers to a Scope which has a common ancestor with prerequisiteScope
        * and that common ancestor is a Scatter block, then find the shard with the same index
        * as 'entry'.  In other words, if you're in the same scatter block as your pre-requisite
        * scope, then depend on the shard (with same index).
        *
        * NOTE: this algorithm was designed for ONE-LEVEL of scattering and probably does not
        * work as-is for nested scatter blocks
        */
      case Some(ancestor: Scatter) => doneKeys.contains(prerequisiteScope.fullyQualifiedName -> entry.index)

      /*
        * Otherwise, simply refer to the collector entry.  This means that 'entry' depends
        * on every shard of the pre-requisite scope to finish.
        */
      case _ => doneKeys.contains(prerequisiteScope.fullyQualifiedName -> None)
    }
  }
}
