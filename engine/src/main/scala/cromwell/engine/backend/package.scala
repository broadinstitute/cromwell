package cromwell.engine

import wdl4s.{Call, Scope}

import scala.concurrent.{ExecutionContext, Future}
import scala.collection.immutable.Seq

package object backend {
  class WorkflowContext(val root: String)
  class CallContext(override val root: String, val stdout: String, val stderr: String) extends WorkflowContext(root)

  final case class ExecutionHash(overallHash: String, dockerHash: Option[String])

  trait JobKey {
    def scope: Scope
    def index: Option[Int]
    def attempt: Int
    def tag: String
  }

  // Ordered by shards, and then ordered by attempts
  type AttemptedCallLogs = Seq[Seq[CallLogs]]
  // Grouped by FQNS
  type WorkflowLogs = Map[FullyQualifiedName, AttemptedCallLogs]

  implicit class EnhancedFutureFuture[A](val ffa: Future[Future[A]])(implicit ec: ExecutionContext) {
    def flatten: Future[A] = ffa flatMap identity
  }

  implicit class EnhancedExecutionHandle(val handle: ExecutionHandle) extends AnyVal {
    def future = Future.successful(handle)
  }

  implicit class EnhancedExecutionResult(val result: ExecutionResult) extends AnyVal {
    def future = Future.successful(result)
  }
}
