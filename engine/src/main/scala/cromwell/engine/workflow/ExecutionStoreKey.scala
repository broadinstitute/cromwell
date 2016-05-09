package cromwell.engine.workflow

import cromwell.backend.JobKey
import cromwell.engine.ExecutionStatus
import cromwell.engine.finalcall.OldStyleFinalCall
import cromwell.engine.workflow.OldStyleWorkflowActor.ExecutionStore
import wdl4s._

import scala.language.postfixOps

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
sealed trait ExecutionStoreKey extends JobKey {
  def tag: String = {
    val shard = index.map(x => s":$x").getOrElse("")
    val attemptTag = if (attempt == 1) "" else s":$attempt"
    s"${scope.unqualifiedName}$shard$attemptTag"
  }
  def retryClone: ExecutionStoreKey
}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
sealed trait OutputKey extends ExecutionStoreKey
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
sealed trait CallKey extends OutputKey
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class BackendCallKey(scope: Call, index: Option[Int], attempt: Int) extends CallKey {
  def retryClone = this.copy(attempt = this.attempt + 1)
}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class CollectorKey(scope: Call, attempt: Int = 1) extends OutputKey {
  override val index: Option[Int] = None
  def retryClone = this.copy(attempt = this.attempt + 1)
}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class ScatterKey(scope: Scatter, index: Option[Int], attempt: Int = 1) extends ExecutionStoreKey {

  /**
   * Creates a sub-ExecutionStore with Starting entries for each of the scoped children.
   *
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
        val shards = (0 until count) map { i => BackendCallKey(call, Option(i), 1) }
        shards :+ CollectorKey(call)
      case scatter: Scatter =>
        throw new UnsupportedOperationException("Nested Scatters are not supported (yet).")
      case e =>
        throw new UnsupportedOperationException(s"Scope ${e.getClass.getName} is not supported.")
    }
  }

  def retryClone = this.copy(attempt = this.attempt + 1)
}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class FinalCallKey(scope: OldStyleFinalCall, attempt: Int = 1) extends CallKey {
  override val index = None
  def retryClone = this.copy(attempt = this.attempt + 1)
}
