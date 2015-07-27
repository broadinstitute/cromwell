package cromwell.engine.backend

import cromwell.binding.values.WdlValue

import scala.concurrent.Future
import scala.util.Try

/**
 * Wraps the execution and cancelling of commands.
 * Creation of the CommandExecution object implies that the execution has already been started.
 */
trait CommandExecution {
  type CallOutputs = Map[String, WdlValue]

  def futureOutputs(): Future[Try[CallOutputs]]

  def cancel(): Boolean
}