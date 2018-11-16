package cromwell.engine.workflow.lifecycle

import akka.actor.FSM

import scala.concurrent.duration._

trait TimedFSM[S] { this: FSM[S, _] =>
  private var lastTransitionTime = System.currentTimeMillis()

  def onTimedTransition(from: S, to: S, duration: FiniteDuration): Unit = {}

  def currentStateDuration: FiniteDuration = (System.currentTimeMillis() - lastTransitionTime).milliseconds

  onTransition {
    case from -> to =>
      val now = System.currentTimeMillis()
      onTimedTransition(from, to, (now - lastTransitionTime).milliseconds)
      lastTransitionTime = now
  }
}
