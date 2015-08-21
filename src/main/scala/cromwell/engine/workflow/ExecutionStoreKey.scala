package cromwell.engine.workflow

import cromwell.binding._
import cromwell.engine.ExecutionStatus
import cromwell.engine.workflow.WorkflowActor.ExecutionStore

sealed trait ExecutionStoreKey {
  val scope: Scope

  val index: Option[Int]

  val parent: Option[ExecutionStoreKey]
}


case class CallKey(override val scope: Call,
                   override val index: Option[Int],
                   override val parent: Option[ExecutionStoreKey])
  extends ExecutionStoreKey


case class ScatterKey(override val scope: Scatter,
                      override val index: Option[Int],
                      //override val status: ExecutionStatus,
                      override val parent: Option[ExecutionStoreKey])
  extends ExecutionStoreKey {

  /**
   * Creates a sub-ExecutionStore with Starting entries for each of the scoped children.
   * @param count Number of ways to scatter the children.
   * @return ExecutionStore of scattered children.
   */
  def populate(count: Int): ExecutionStore = {
    val keys = for {
      childScope <- this.scope.children
      collector = collectorForChild(childScope, count)
      entry <- collector +: collector.keys
    } yield entry

    (keys map {
      _ -> ExecutionStatus.Starting
    }).toMap
  }

  private def collectorForChild(scope: Scope, count: Int): CollectorKey = {
    val parent = Option(this)
    val indexed = (0 until count) map { i =>
      val index = Option(i)
      scope match {
        case call: Call =>
          CallKey(call, index, parent)
        case scatter: Scatter =>
          ScatterKey(scatter, index, parent)
      }
    }
    CollectorKey(scope, parent, indexed)
  }
}


case class CollectorKey(override val scope: Scope,
                        override val parent: Option[ExecutionStoreKey],
                        keys: Seq[ExecutionStoreKey])
  extends ExecutionStoreKey {
  override val index = None
}
