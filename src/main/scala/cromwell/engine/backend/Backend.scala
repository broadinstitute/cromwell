package cromwell.engine.backend

import cromwell.binding.WdlExpression.ScopedLookupFunction
import cromwell.binding.values.WdlValue
import cromwell.binding.{Call, TaskOutput}

import scala.util.Try

/**
 * trait to be implemented by concrete backends.
 */
trait Backend {

  /**
   * Execute the specified command line using the provided symbol store, evaluating the task outputs to produce
   * a mapping of local task output names to WDL values.
   */
  def executeCommand(commandLine: String, call: Call, taskOutputs: Seq[TaskOutput], scopedLookupFunction: ScopedLookupFunction): Map[String, Try[WdlValue]]

}
