package cromwell.engine.workflow.lifecycle

package execution {

  import cromwell.core.CallKey
  import wom.values.WomEvaluatedCallInputs

  final case class JobRunning(key: CallKey, inputs: WomEvaluatedCallInputs)
  final case class JobStarting(callKey: CallKey)
}
