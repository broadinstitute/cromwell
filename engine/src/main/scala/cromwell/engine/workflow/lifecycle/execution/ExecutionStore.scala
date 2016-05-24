package cromwell.engine.workflow.lifecycle.execution

import akka.actor.ActorRef
import cromwell.backend.{BackendJobDescriptorKey, JobKey}
import cromwell.engine.ExecutionStatus
import cromwell.engine.ExecutionStatus._
import cromwell.engine.workflow.lifecycle.execution.ExecutionStore.ExecutionStoreEntry
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.CollectorKey
import wdl4s.{Scatter, Scope}

import scala.language.postfixOps

object ExecutionStore {
  def empty = ExecutionStore(Map.empty)
  type ExecutionStoreEntry = (JobKey, ExecutionStatus)
}

case class ExecutionStore(store: Map[JobKey, ExecutionStatus]) {
  def add(values: Map[JobKey, ExecutionStatus]) = {
    values foreach { case (key, value) => System.out.println(s"Adding $key -> $value") }
    this.copy(store = store ++ values)
  }

  /** Find currently runnable scopes */
  def runnableScopes: Iterable[JobKey] = store.filter(isRunnable).keys

  private def isRunnable(entry: ExecutionStoreEntry) = {
    entry match {
      case (key, ExecutionStatus.NotStarted) => arePrerequisitesDone(key)
      case _ => false
    }
  }

  def findShardEntries(key: CollectorKey): Iterable[ExecutionStoreEntry] = store collect {
    case (k: BackendJobDescriptorKey, v) if k.scope == key.scope && k.isShard => (k, v)
  }

  private def arePrerequisitesDone(key: JobKey): Boolean = {
    def isDone(e: JobKey): Boolean = store exists { case (k, s) => k.scope == e.scope && k.index == e.index && s == ExecutionStatus.Done }

    val upstream = key.scope.prerequisiteScopes.map(s => upstreamEntries(key, s))
    val downstream = key match {
      case collector: CollectorKey => findShardEntries(collector)
      case _ => Nil
    }
    val dependencies = upstream.flatten ++ downstream
    val dependenciesResolved = store.filter(dependencies).keys forall isDone

    /**
      * We need to make sure that all prerequisiteScopes have been resolved to some entry before going forward.
      * If a scope cannot be resolved it may be because it is in a scatter that has not been populated yet,
      * therefore there is no entry in the executionStore for this scope.
      * If that's the case this prerequisiteScope has not been run yet, hence the (upstream forall {_.nonEmpty})
      */
    (upstream forall { _.nonEmpty }) && dependenciesResolved
  }

  private def upstreamEntries(entry: JobKey, prerequisiteScope: Scope): Seq[ExecutionStoreEntry] = {
    prerequisiteScope.closestCommonAncestor(entry.scope) match {
      /**
        * If this entry refers to a Scope which has a common ancestor with prerequisiteScope
        * and that common ancestor is a Scatter block, then find the shard with the same index
        * as 'entry'.  In other words, if you're in the same scatter block as your pre-requisite
        * scope, then depend on the shard (with same index).
        *
        * NOTE: this algorithm was designed for ONE-LEVEL of scattering and probably does not
        * work as-is for nested scatter blocks
        */
      case Some(ancestor: Scatter) =>
        store filter { case (k, _) => k.scope == prerequisiteScope && k.index == entry.index } toSeq

      /**
        * Otherwise, simply refer to the entry the collector entry.  This means that 'entry' depends
        * on every shard of the pre-requisite scope to finish.
        */
      case _ =>
        store filter { case (k, _) => k.scope == prerequisiteScope && k.index.isEmpty } toSeq
    }
  }
}
