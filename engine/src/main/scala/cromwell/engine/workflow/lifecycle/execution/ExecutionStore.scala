package cromwell.engine.workflow.lifecycle.execution

import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.ExecutionStatus._
import cromwell.core.{CallKey, ExecutionStatus, JobKey}
import cromwell.engine.workflow.lifecycle.execution.ExecutionStore.ExecutionStoreEntry
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.{CollectorKey, DeclarationKey, ScatterKey, SubWorkflowKey}
import wdl4s._


object ExecutionStore {
  def empty = ExecutionStore(Map.empty[JobKey, ExecutionStatus])
  type ExecutionStoreEntry = (JobKey, ExecutionStatus)
  def apply(workflow: Workflow) = {
    // Only add direct children to the store, the rest is dynamically created when necessary
    val keys = workflow.children map {
      case call: TaskCall => Option(BackendJobDescriptorKey(call, None, 1))
      case call: WorkflowCall => Option(SubWorkflowKey(call, None, 1))
      case scatter: Scatter => Option(ScatterKey(scatter))
      case declaration: Declaration if declaration.expression.isDefined => Option(DeclarationKey(declaration, None))
      case _ => None // Ifs will need to be added here when supported
    }

    new ExecutionStore(keys.flatten.map(_ -> NotStarted).toMap)
  }
}

case class ExecutionStore(store: Map[JobKey, ExecutionStatus]) {
  def add(values: Map[JobKey, ExecutionStatus]) = this.copy(store = store ++ values)

  // Convert the store to a `List` before `collect`ing to sidestep expensive and pointless hashing of `Scope`s when
  // assembling the result.
  def runnableScopes = store.toList collect { case entry if isRunnable(entry) => entry._1 }

  def backendJobKeys = store.keys.toList collect { case k: BackendJobDescriptorKey => k }

  private def isRunnable(entry: ExecutionStoreEntry) = {
    entry match {
      case (key, ExecutionStatus.NotStarted) => arePrerequisitesDone(key)
      case _ => false
    }
  }

  def findShardEntries(key: CollectorKey): List[ExecutionStoreEntry] = store.toList filter {
    case (k: CallKey, v) => k.scope == key.scope && k.isShard
    case _ => false
  }

  private def arePrerequisitesDone(key: JobKey): Boolean = {
    val upstream = key.scope match {
      case node: GraphNode => node.upstream collect {
        case n: Call => upstreamEntry(key, n)
        case n: Scatter => upstreamEntry(key, n)
        case n: Declaration if n.expression.isDefined => upstreamEntry(key, n)
      }
      case _ => Set.empty
    }

    val downstream: List[(JobKey, ExecutionStatus)] = key match {
      case collector: CollectorKey => findShardEntries(collector)
      case _ => Nil
    }

    val dependencies = upstream.flatten ++ downstream
    val dependenciesResolved = dependencies forall { case (_, s) => s == ExecutionStatus.Done }

    /*
      * We need to make sure that all prerequisiteScopes have been resolved to some entry before going forward.
      * If a scope cannot be resolved it may be because it is in a scatter that has not been populated yet,
      * therefore there is no entry in the executionStore for this scope.
      * If that's the case this prerequisiteScope has not been run yet, hence the (upstream forall {_.nonEmpty})
      */
    (upstream forall { _.nonEmpty }) && dependenciesResolved
  }

  private def upstreamEntry(entry: JobKey, prerequisiteScope: Scope): Option[ExecutionStoreEntry] = {
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
        store find {
          case (k, _) => k.scope == prerequisiteScope && k.index == entry.index
        }

      /*
        * Otherwise, simply refer to the collector entry.  This means that 'entry' depends
        * on every shard of the pre-requisite scope to finish.
        */
      case _ =>
        store find {
          case (k, _) => k.scope == prerequisiteScope && k.index.isEmpty
        }
    }
  }
}
