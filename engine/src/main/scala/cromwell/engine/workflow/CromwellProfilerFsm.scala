package cromwell.engine.workflow

import akka.actor.LoggingFSM

/* Super class of the messages that are published to the Akka event stream */
trait CromwellFsmStateAndData

/**
  * Super class of the FSM actors implemented in Cromwell that wish to publish their FSM states and data.
  * SubClasses that extend from this class need to implement the `wrappedStateAndData` method, and return an ADT
  * that is a subtype of `CromwellFsmStateAndData`.
  * This event is then published in the Akka's actor system event stream.
  *
  * @tparam S StateName to be used with the FSM
  * @tparam D StateData to be used with the FSM
  */
trait CromwellProfilerFsm[S,D] extends LoggingFSM[S,D] {

  def wrappedStateAndData[? <: CromwellFsmStateAndData](state: S, data: D): CromwellFsmStateAndData

  onTransition {
    case _ -> newState =>
      context.system.eventStream.publish(wrappedStateAndData(newState, nextStateData))
  }

}