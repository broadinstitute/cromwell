package cromwell.engine.workflow

import cromwell.binding._
import cromwell.engine.ExecutionStatus
import cromwell.engine.workflow.WorkflowActor.ExecutionStore

sealed trait ExecutionStoreKey {
  val scope: Scope

  val index: Option[Int]

  val parent: Option[ExecutionStoreKey]
}

trait OutputKey extends ExecutionStoreKey

case class CallKey(scope: Call, index: Option[Int], parent: Option[ExecutionStoreKey]) extends OutputKey
case class CollectorKey(scope: Call, parent: Option[ExecutionStoreKey]) extends OutputKey {
  val index: Option[Int] = None
}

case class ScatterKey(scope: Scatter, index: Option[Int], parent: Option[ExecutionStoreKey]) extends ExecutionStoreKey {

  /**
   * Creates a sub-ExecutionStore with Starting entries for each of the scoped children.
   * @param count Number of ways to scatter the children.
   * @return ExecutionStore of scattered children.
   */
  def populate(count: Int): ExecutionStore = {
    val keys = this.scope.children flatMap { explode(_, count) }

    (keys map {
      _ -> ExecutionStatus.NotStarted
    }).toMap
  }

  private def explode(scope: Scope, count: Int): Seq[ExecutionStoreKey] = {
    val parent = Option(this)
      scope match {
        case call: Call =>
          val shards = (0 until count) map { i => CallKey(call, Option(i), parent) }
          shards :+ CollectorKey(call, parent)
        case scatter: Scatter =>
          throw new UnsupportedOperationException("Nested Scatters are not supported (yet).")
        case e =>
          throw new UnsupportedOperationException(s"Scope ${e.getClass.getName} is not supported.")
    }
  }
}

