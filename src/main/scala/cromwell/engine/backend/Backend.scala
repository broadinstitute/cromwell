package cromwell.engine.backend

import com.typesafe.config.Config
import cromwell.binding._
import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine._
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.backend.jes.JesBackend
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.backend.sge.SgeBackend
import cromwell.engine.db.DataAccess
import cromwell.engine.workflow.CallKey
import cromwell.parser.BackendType

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

object Backend {
  class StdoutStderrException(message: String) extends RuntimeException(message)
  def from(backendConf: Config): Backend = from(backendConf.getString("backend"))
  def from(name: String) = name.toLowerCase match {
    case "local" => new LocalBackend
    case "jes" => new JesBackend
    case "sge" => new SgeBackend
    case doh => throw new IllegalArgumentException(s"$doh is not a recognized backend")
  }
  case class RestartableWorkflow(id: WorkflowId, source: WorkflowSourceFiles)
}

/**
 * Trait to be implemented by concrete backends.
 */
trait Backend {
  type BackendCall <: backend.BackendCall

  /**
   * Return a possibly altered copy of inputs reflecting any localization of input file paths that might have
   * been performed for this `Backend` implementation.
   */
  def adjustInputPaths(call: Call, inputs: CallInputs): CallInputs

  def adjustOutputPaths(call: Call, outputs: CallOutputs): CallOutputs

  /**
   * Do whatever work is required to initialize the workflow, returning a copy of
   * the coerced inputs present in the `WorkflowDescriptor` with any input `WdlFile`s
   * adjusted for the host workflow execution path.
   */
  def initializeForWorkflow(workflow: WorkflowDescriptor, dataAccess: DataAccess): Try[HostInputs]

  /**
   * Do whatever cleaning up work is required when a workflow reaches a terminal state.
   */
  def cleanUpForWorkflow(workflow: WorkflowDescriptor)(implicit ec: ExecutionContext): Future[Any]

  /**
   * Execute the Call (wrapped in a BackendCall), return the outputs if it is
   * successful, otherwise, returns Failure with a reason why the execution failed
   */
  def execute(backendCall: BackendCall): ExecutionResult

  /**
   * Essentially turns a Call object + CallInputs into a BackendCall
   */
  def bindCall(workflowDescriptor: WorkflowDescriptor,
               key: CallKey,
               locallyQualifiedInputs: CallInputs,
               abortRegistrationFunction: AbortRegistrationFunction): BackendCall

  /**
   * Do whatever is appropriate for this backend implementation to support restarting the specified workflows.
   */
  def handleCallRestarts(restartableWorkflows: Seq[RestartableWorkflow], dataAccess: DataAccess)(implicit ec: ExecutionContext): Future[Any]

  /**
   * Return CallStandardOutput which contains the stdout/stderr of the particular call
   */
  def stdoutStderr(descriptor: WorkflowDescriptor, callName: String, index: ExecutionIndex): StdoutStderr

  def backendType: BackendType

  /**
   * Validate that workflow options contain all required information.
   */
  @throws[IllegalArgumentException]("if a value is missing / incorrect")
  def assertWorkflowOptions(options: Map[String, String])

  def makeTag(backendCall: BackendCall): String =
    s"${this.getClass.getSimpleName} [UUID(${backendCall.workflowDescriptor.shortId}):${backendCall.call.name}]"
}
