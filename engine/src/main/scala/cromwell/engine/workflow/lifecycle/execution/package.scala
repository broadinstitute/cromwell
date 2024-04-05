package cromwell.engine.workflow.lifecycle

package execution {

  import cromwell.core.CallKey
  import wom.values.WomEvaluatedCallInputs

  import scala.util.control.NoStackTrace

  final case class JobRunning(key: CallKey, inputs: WomEvaluatedCallInputs)
  final case class JobStarting(callKey: CallKey)

  /**
    * An exception specific to conditions inside the executing WDL, as opposed to one that is "Cromwell's fault"
    * @param message Description suitable for user display
    */
  final case class WdlRuntimeException(message: String) extends RuntimeException with NoStackTrace {
    override def getMessage: String = message
  }
}
