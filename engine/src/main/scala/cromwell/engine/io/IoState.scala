package cromwell.engine.io

import cats.Show
import cromwell.core.CromwellFatalExceptionMarker
import cromwell.core.retry.IORetry.StatefulIoError
import org.apache.commons.lang3.exception.ExceptionUtils

object IoState {
  object EnhancedCromwellIoException {
    def apply[S](state: S, cause: Throwable)(implicit showState: Show[S]): EnhancedCromwellIoException = {
      EnhancedCromwellIoException(s"[${showState.show(state)}] - ${ExceptionUtils.getMessage(cause)}", cause)
    }
  }
  
  case class EnhancedCromwellIoException(message: String, cause: Throwable) 
    extends Throwable(message, cause, true, false) with CromwellFatalExceptionMarker
  
  implicit val showState = new Show[IoState] {
    override def show(t: IoState) = s"Attempted ${t.attempt} time(s)"
  }
  
  implicit val stateToThrowable = new StatefulIoError[IoState] {
    override def toThrowable(state: IoState, throwable: Throwable) = {
      state.throwables.foreach(throwable.addSuppressed)
      EnhancedCromwellIoException(state, throwable)
    }
  }

  val updateState: (Throwable, IoState) => IoState = (throwable, state) => {
    state.copy(attempt = state.attempt + 1, throwables = state.throwables :+ throwable)
  }
}

case class IoState(attempt: Int, throwables: List[Throwable] = List.empty)
