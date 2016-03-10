package cromwell.engine

import scala.concurrent.{ExecutionContext, Future}
import scala.collection.immutable.Seq

package object backend {
  class WorkflowContext(val root: String)
  class CallContext(override val root: String, val stdout: String, val stderr: String) extends WorkflowContext(root)

  // Ordered by shards, and then ordered by attempts
  type AttemptedCallLogs = Seq[Seq[CallLogs]]
  // Grouped by FQNS
  type WorkflowLogs = Map[FullyQualifiedName, AttemptedCallLogs]

  implicit class EnhancedFutureFuture[A](val ffa: Future[Future[A]])(implicit ec: ExecutionContext) {
    def flatten: Future[A] = ffa flatMap { fa => fa }
  }

  implicit class EnhancedExecutionHandle(val handle: ExecutionHandle) extends AnyVal {
    def future = Future.successful(handle)
  }

  implicit class EnhancedExecutionResult(val result: ExecutionResult) extends AnyVal {
    def future = Future.successful(result)
  }
}
