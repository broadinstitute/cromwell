package cromwell.engine.io

import cats.Show

object IoState {
  implicit val showState = new Show[IoState] {
    override def show(state: IoState) = s"Attempted ${state.attempt} time(s)"
  }

  val updateState: (Throwable, IoState) => IoState = (_, state) => {
    state.copy(attempt = state.attempt + 1)
  }
}

case class IoState(attempt: Int)
