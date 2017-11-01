package cromwell.engine.workflow.lifecycle

import akka.actor.ActorRef

package execution {

  import cromwell.core.CallKey
  import wom.values.WomEvaluatedCallInputs

  final case class JobRunning(key: CallKey, inputs: WomEvaluatedCallInputs, executionActor: Option[ActorRef])
  final case class JobStarting(callKey: CallKey)
}
