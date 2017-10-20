package cromwell.backend

import scala.concurrent.{ExecutionContext, Future}


package object async {
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
