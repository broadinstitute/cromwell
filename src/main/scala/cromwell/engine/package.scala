package cromwell

import cromwell.binding._
import org.joda.time.DateTime

import scala.language.implicitConversions
import scalaz.ValidationNel

package object engine {

  class WorkflowContext(val root: String)
  class CallContext(override val root: String, val stdout: String, val stderr: String) extends WorkflowContext(root)

  /**
   * Represents the collection of source files that a user submits to run a workflow
   */
  final case class WorkflowSourceFiles(wdlSource: WdlSource, inputsJson: WdlJson, workflowOptionsJson: WorkflowOptionsJson)

  final case class AbortFunction(function: ()=>Unit)
  final case class AbortRegistrationFunction(register: AbortFunction=>Unit)

  final case class ExecutionEventEntry(description: String, startTime: DateTime, endTime: DateTime)
  final case class ExecutionHash(overallHash: String, dockerHash: Option[String])

  type ErrorOr[+A] = ValidationNel[String, A]
}
