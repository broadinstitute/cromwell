package centaur.reporting

import cats.effect.IO
import centaur.test.submit.{SubmitResponse, SubmitWorkflowResponse}

trait SuccessReporter {
  def logSuccessfulRun(submitResponse: SubmitWorkflowResponse): IO[Unit]
}

object SuccessReporters {
  /*
   * This is gross and piggy backs on the error reporting code but achieves what we need for now without a refactoring
   * of the error reporting code to handle both success and error reporting
   */
  private val successReporters: List[ErrorReporter with SuccessReporter] =
    ErrorReporters.errorReporters.errorReporters.collect { case s: SuccessReporter => s }

  def logSuccessfulRun(submitResponse: SubmitWorkflowResponse): IO[SubmitResponse] =
    if (successReporters.isEmpty) IO.pure(submitResponse)
    else {
      val listIo = successReporters.map(_.logSuccessfulRun(submitResponse))
      AggregatedIo.aggregateExceptions("Errors while reporting centaur success", listIo).map(_ => submitResponse)
    }
}
