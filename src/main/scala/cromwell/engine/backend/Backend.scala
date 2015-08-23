package cromwell.engine.backend

import java.nio.file.{Files, Paths}

import com.typesafe.config.Config
import cromwell.binding
import cromwell.binding.WdlExpression.ScopedLookupFunction
import cromwell.binding._
import cromwell.binding.types.{WdlFileType, WdlMapType}
import cromwell.binding.values.{WdlFile, WdlMap, WdlString, WdlValue}
import cromwell.engine._
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.backend.jes.JesBackend
import cromwell.engine.backend.local.{LocalBackendCall, LocalBackend, LocalEngineFunctions}
import cromwell.engine.backend.sge.SgeBackend
import cromwell.engine.db.DataAccess
import cromwell.parser.BackendType

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

object Backend {
  class StdoutStderrException(message: String) extends RuntimeException(message)
  def from(backendConf: Config): Backend = {
    backendConf.getString("backend").toLowerCase match {
      case "local" => new LocalBackend
      case "jes" => new JesBackend
      case "sge" => new SgeBackend
      case doh => throw new IllegalArgumentException(s"$doh is not a recognized backend")
    }
  }

  case class RestartableWorkflow(id: WorkflowId, source: WdlSource, json: WdlJson, inputs: binding.WorkflowRawInputs)
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
  def initializeForWorkflow(workflow: WorkflowDescriptor): Try[HostInputs]

  /**
   * Execute the Call (wrapped in a BackendCall), return the outputs if it is
   * successful, otherwise, returns Failure with a reason why the execution failed
   */
  def execute(backendCall: BackendCall): Try[CallOutputs]

  /**
   * Essentially turns a Call object + CallInputs into a BackendCall
   */
  def bindCall(workflowDescriptor: WorkflowDescriptor,
               call: Call,
               locallyQualifiedInputs: CallInputs,
               abortRegistrationFunction: AbortRegistrationFunction): BackendCall

  /**
   * Do whatever is appropriate for this backend implementation to support restarting the specified workflows.
   */
  def handleCallRestarts(restartableWorkflows: Seq[RestartableWorkflow], dataAccess: DataAccess)(implicit ec: ExecutionContext): Future[Any]

  /**
   * Return CallStandardOutput which contains the stdout/stderr of the particular call
   */
  def stdoutStderr(workflowId: WorkflowId, workflowName: String, callName: String): StdoutStderr

  def backendType: BackendType

  def makeTag(backendCall: BackendCall): String =
    s"${this.getClass.getSimpleName} [UUID(${backendCall.workflowDescriptor.shortId}):${backendCall.call.name}]"
}
