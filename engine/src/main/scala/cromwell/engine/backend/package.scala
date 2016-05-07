package cromwell.engine

import scala.collection.immutable.Seq
import scala.concurrent.{ExecutionContext, Future}

package object backend {
  // Ordered by shards, and then ordered by attempts
  type AttemptedCallLogs = Seq[Seq[CallLogs]]
  // Grouped by FQNS
  type WorkflowLogs = Map[FullyQualifiedName, AttemptedCallLogs]

  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  implicit class EnhancedFutureFuture[A](val ffa: Future[Future[A]])(implicit ec: ExecutionContext) {
    def flatten: Future[A] = ffa flatMap identity
  }

  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  implicit class EnhancedExecutionHandle(val handle: OldStyleExecutionHandle) extends AnyVal {
    def future = Future.successful(handle)
  }

  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  implicit class EnhancedExecutionResult(val result: OldStyleExecutionResult) extends AnyVal {
    def future = Future.successful(result)
  }
}
