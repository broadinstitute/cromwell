package cromwell.engine

import wdl4s.{Call, Scope}

import scala.concurrent.{ExecutionContext, Future}
import scala.collection.immutable.Seq

package object backend {
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  final case class ExecutionHash(overallHash: String, dockerHash: Option[String])

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
