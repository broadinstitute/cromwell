package cromwell.engine.workflow

import cromwell.engine.finalcall.FinalCall
import wdl4s._
import cromwell.engine.{WorkflowDescriptor, ExecutionStatus}
import cromwell.engine.workflow.WorkflowActor.ExecutionStore

import scala.language.postfixOps

sealed trait ExecutionStoreKey {
  def scope: Scope
  def index: Option[Int]
  def tag: String = {
    val shard = index.map(x => s":$x").getOrElse("")
    val attemptTag = if (attempt == 1) "" else s":$attempt"
    s"${scope.unqualifiedName}$shard$attemptTag"
  }
  def attempt: Int
  def retryClone: ExecutionStoreKey
}

sealed trait OutputKey extends ExecutionStoreKey

sealed trait CallKey extends OutputKey

object BackendCallKey{
  private val CallPrefix = "call"
  private val ShardPrefix = "shard"
  private val AttemptPrefix = "attempt"
}
case class BackendCallKey(scope: Call, index: Option[Int], attempt: Int = 1) extends CallKey {
  import BackendCallKey._

  def retryClone = this.copy(attempt = this.attempt + 1)

  def callRootPathWithBaseRoot(descriptor: WorkflowDescriptor, baseRoot: String) = {
    val call = s"$CallPrefix-${scope.unqualifiedName}"
    val shard = index map { s => s"$ShardPrefix-$s" } getOrElse ""
    val retry = if (attempt > 1) s"$AttemptPrefix-$attempt" else ""
    descriptor.workflowRootPathWithBaseRoot(baseRoot).resolve(call).resolve(shard).resolve(retry)
  }
}
case class CollectorKey(scope: Call, attempt: Int = 1) extends OutputKey {
  override val index: Option[Int] = None
  def retryClone = this.copy(attempt = this.attempt + 1)
}

case class ScatterKey(scope: Scatter, index: Option[Int], attempt: Int = 1) extends ExecutionStoreKey {

  /**
   * Creates a sub-ExecutionStore with Starting entries for each of the scoped children.
   * @param count Number of ways to scatter the children.
   * @return ExecutionStore of scattered children.
   */
  def populate(count: Int): ExecutionStore = {
    val keys = this.scope.children flatMap { explode(_, count) }
    keys map { _ -> ExecutionStatus.NotStarted } toMap
  }

  private def explode(scope: Scope, count: Int): Seq[ExecutionStoreKey] = {
    scope match {
      case call: Call =>
        val shards = (0 until count) map { i => BackendCallKey(call, Option(i)) }
        shards :+ CollectorKey(call)
      case scatter: Scatter =>
        throw new UnsupportedOperationException("Nested Scatters are not supported (yet).")
      case e =>
        throw new UnsupportedOperationException(s"Scope ${e.getClass.getName} is not supported.")
    }
  }

  def retryClone = this.copy(attempt = this.attempt + 1)
}

case class FinalCallKey(scope: FinalCall, attempt: Int = 1) extends CallKey {
  override val index = None
  def retryClone = this.copy(attempt = this.attempt + 1)
}
