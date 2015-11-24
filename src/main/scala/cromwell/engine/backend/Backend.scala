package cromwell.engine.backend

import com.typesafe.config.Config
import cromwell.binding._
import cromwell.binding.expression.WdlStandardLibraryFunctions
import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine._
import cromwell.engine.backend.jes.JesBackend
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.backend.sge.SgeBackend
import cromwell.engine.db.ExecutionDatabaseKey
import cromwell.engine.workflow.{CallKey, WorkflowOptions}
import cromwell.logging.WorkflowLogger
import cromwell.parser.BackendType
import org.slf4j.LoggerFactory

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

object Backend {
  lazy val LocalBackend = new LocalBackend
  lazy val JesBackend = new JesBackend { jesConf } //forces configuration resolution to fail now if something is missing
  lazy val SgeBackend = new SgeBackend

  class StdoutStderrException(message: String) extends RuntimeException(message)

  def from(backendConf: Config): Backend = from(backendConf.getString("backend"))
  def from(name: String) = name.toLowerCase match {
    case "local" => LocalBackend
    case "jes" => JesBackend
    case "sge" => SgeBackend
    case doh => throw new IllegalArgumentException(s"$doh is not a recognized backend")
  }

  case class RestartableWorkflow(id: WorkflowId, source: WorkflowSourceFiles)
}

trait JobKey

/**
 * Trait to be implemented by concrete backends.
 */
trait Backend {

  type BackendCall <: backend.BackendCall
  type IOInterface <: cromwell.binding.IOInterface

  /**
   * Return a possibly altered copy of inputs reflecting any localization of input file paths that might have
   * been performed for this `Backend` implementation.
   */
  def adjustInputPaths(callKey: CallKey, inputs: CallInputs, workflowDescriptor: WorkflowDescriptor): CallInputs

  // FIXME: This is never called...
  def adjustOutputPaths(call: Call, outputs: CallOutputs): CallOutputs

  /**
   * Do whatever work is required to initialize the workflow, returning a copy of
   * the coerced inputs present in the `WorkflowDescriptor` with any input `WdlFile`s
   * adjusted for the host workflow execution path.
   */
  def initializeForWorkflow(workflow: WorkflowDescriptor): Try[HostInputs]

  /**
   * Do whatever cleaning up work is required when a workflow reaches a terminal state.
   */
  def cleanUpForWorkflow(workflow: WorkflowDescriptor)(implicit ec: ExecutionContext): Future[Any] = Future.successful({})

  /**
   * Essentially turns a Call object + CallInputs into a BackendCall
   */
  def bindCall(workflowDescriptor: WorkflowDescriptor,
               key: CallKey,
               locallyQualifiedInputs: CallInputs,
               abortRegistrationFunction: AbortRegistrationFunction): BackendCall

  /**
   * Engine functions that don't need a Call context (e.g. read_lines(), read_float(), etc)
   */
  def engineFunctions(interface: IOInterface): WdlStandardLibraryFunctions

  /**
    * Interface to be used primarily by engine functions requiring IO capabilities.
    */
  def ioInterface(workflowOptions: WorkflowOptions): IOInterface

  /**
   * Do any work that needs to be done <b>before</b> attempting to restart a workflow.
   */
  def prepareForRestart(restartableWorkflow: WorkflowDescriptor)(implicit ec: ExecutionContext): Future[Unit]

  /**
   * Return CallStandardOutput which contains the stdout/stderr of the particular call
   */
  def stdoutStderr(descriptor: WorkflowDescriptor, callName: String, index: ExecutionIndex): CallLogs

  /**
   * Provides a function that given a WdlFile, returns its hash.
   */
  def fileHasher(workflowDescriptor: WorkflowDescriptor): FileHasher

  def backendType: BackendType

  /**
   * Validate that workflow options contain all required information.
   */
  @throws[IllegalArgumentException]("if a value is missing / incorrect")
  def assertWorkflowOptions(options: WorkflowOptions): Unit = {}

  private def backendClass = backendType.toString.toLowerCase.capitalize + "Backend"

  /** Default implementation assumes backends do not support resume, returns an empty Map. */
  def findResumableExecutions(id: WorkflowId)(implicit ec: ExecutionContext): Future[Map[ExecutionDatabaseKey, JobKey]] = Future.successful(Map.empty)

  def workflowLogger(descriptor: WorkflowDescriptor) = WorkflowLogger(
    backendClass,
    descriptor,
    otherLoggers = Seq(LoggerFactory.getLogger(getClass.getName))
  )

  def workflowLoggerWithCall(backendCall: BackendCall) = WorkflowLogger(
    backendClass,
    backendCall.workflowDescriptor,
    otherLoggers = Seq(LoggerFactory.getLogger(getClass.getName)),
    callTag = Option(backendCall.key.tag)
  )
}
