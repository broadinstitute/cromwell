package cromwell.engine.workflow.lifecycle.execution

import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.ExecutionStatus._
import cromwell.core.{CallKey, JobKey}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.{apply => _, _}
import wdl4s._


object ExecutionStore {
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
}

case class ExecutionStore(private val statusStore: Map[JobKey, ExecutionStatus], hasNewRunnables: Boolean) {

  // View of the statusStore more suited for lookup based on status
  lazy val store: Map[ExecutionStatus, List[JobKey]] = statusStore.groupBy(_._2).mapValues(_.keys.toList)

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

    keysWithStatus(NotStarted).exists(jobKey => !upstreamFailed(jobKey.scope)) ||
      keysWithStatus(QueuedInCromwell).nonEmpty ||
      keysWithStatus(Starting).nonEmpty ||
      keysWithStatus(Running).nonEmpty
  }

  def jobStatus(jobKey: JobKey): Option[ExecutionStatus] = statusStore.get(jobKey)

  def startedJobs: List[BackendJobDescriptorKey] = {
    store.filterNot(_ == NotStarted).values.toList.flatten collect {
      case k: BackendJobDescriptorKey => k
    }
  }

  private def hasFailedScope(s: GraphNode): Boolean = keysWithStatus(Failed).exists(_.scope == s)

  def hasFailedJob: Boolean = keysWithStatus(Failed).nonEmpty

  override def toString = store.map { case (j, s) => s"$j -> $s" } mkString System.lineSeparator()

  def add(values: Map[JobKey, ExecutionStatus]) = {
    val hasTerminalStatus = values.values.exists(_.isTerminal)
    this.copy(statusStore = statusStore ++ values, hasNewRunnables = hasNewRunnables || hasTerminalStatus)
  }

  // Convert the store to a `List` before `collect`ing to sidestep expensive and pointless hashing of `Scope`s when
  // assembling the result.
  def runnableScopes = {
    val doneKeys = store.filterKeys(_.isDoneOrBypassed).values.toList.flatten
    keysWithStatus(NotStarted) filter arePrerequisitesDone(doneKeys)
  }

  def findCompletedShardsForOutput(key: CollectorKey): List[JobKey] = {
    store.filterKeys(_.isDoneOrBypassed).values.toList.flatten collect {
      case k: CallKey if k.scope == key.scope && k.isShard => k
      case k: DynamicDeclarationKey if k.scope == key.scope && k.isShard => k
    }
  }

  // Just used to decide whether a collector can be run. In case the shard entries haven't been populated into the
  // execution store yet.
  private case class TempJobKey(scope: Scope with GraphNode, index: Option[Int]) extends JobKey {
    // If these are ever used, we've done something wrong...
    override def attempt: Int = throw new NotImplementedError("We've done something wrong.")
    override def tag: String = throw new NotImplementedError("We've done something wrong.")
  }
  
  private def emulateShardEntries(key: CollectorKey): List[JobKey] = {
    (0 until key.scatterWidth).toList flatMap { i => key.scope match {
      case c: Call => List(TempJobKey(c, Option(i)))
      case d: Declaration => List(TempJobKey(d, Option(i)))
      case _ => throw new RuntimeException("Don't collect that.")
    }}
  }

  private def arePrerequisitesDone(doneKeys: List[JobKey])(key: JobKey): Boolean = {
    val upstreamAreDone = key.scope.upstream forall {
      case n: Call => upstreamIsDone(key, n, doneKeys)
      case n: Scatter => upstreamIsDone(key, n, doneKeys)
      case n: Declaration => upstreamIsDone(key, n, doneKeys)
      case _ => true
    }

    val shardEntriesForCollectorAreDone: Boolean = key match {
      case collector: CollectorKey => emulateShardEntries(collector).forall(shardKey => doneKeys.exists({ doneKey =>
        shardKey.scope.fullyQualifiedName == doneKey.scope.fullyQualifiedName &&
          shardKey.index == doneKey.index
      }))
      case _ => true
    }

    upstreamAreDone && shardEntriesForCollectorAreDone
  }

  private def upstreamIsDone(entry: JobKey, prerequisiteScope: Scope, doneKeys: List[JobKey]): Boolean = {
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
      case Some(ancestor: Scatter) =>
        doneKeys exists { k =>
          k.scope == prerequisiteScope && k.index == entry.index
        }

      /*
        * Otherwise, simply refer to the collector entry.  This means that 'entry' depends
        * on every shard of the pre-requisite scope to finish.
        */
      case _ =>
        doneKeys exists { k =>
          k.scope == prerequisiteScope && k.index.isEmpty
        }
    }
  }
}
