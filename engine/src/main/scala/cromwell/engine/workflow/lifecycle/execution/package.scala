package cromwell.engine.workflow.lifecycle

import akka.actor.ActorRef
import wdl4s._

package execution {

  import cromwell.core.CallKey

  final case class JobRunning(key: CallKey, inputs: EvaluatedTaskInputs, executionActor: Option[ActorRef])
  final case class JobStarting(callKey: CallKey)
}
