package cromwell.engine.io

import cats.Show
import common.util.IORetry.StatefulIoError
import cromwell.core.CromwellFatalExceptionMarker
import org.apache.commons.lang3.exception.ExceptionUtils

object IoAttempts {
  object EnhancedCromwellIoException {
    def apply[S](state: S, cause: Throwable)(implicit showState: Show[S]): EnhancedCromwellIoException =
      EnhancedCromwellIoException(s"[${showState.show(state)}] - ${ExceptionUtils.getMessage(cause)}", cause)
  }

  case class EnhancedCromwellIoException(message: String, cause: Throwable)
      extends Throwable(message, cause, true, false)
      with CromwellFatalExceptionMarker

  implicit val showState = new Show[IoAttempts] {
    override def show(t: IoAttempts) = s"Attempted ${t.attempts} time(s)"
  }

  implicit val stateToThrowable = new StatefulIoError[IoAttempts] {
    override def toThrowable(state: IoAttempts, throwable: Throwable) = {
      state.throwables.foreach(throwable.addSuppressed)
      EnhancedCromwellIoException(state, throwable)
    }
  }

  val updateState: (Throwable, IoAttempts) => IoAttempts = (throwable, state) =>
    state.copy(attempts = state.attempts + 1, throwables = state.throwables :+ throwable)
}

case class IoAttempts(attempts: Int, throwables: List[Throwable] = List.empty)
