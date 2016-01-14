package cromwell.engine.backend

import akka.actor.ActorSystem
import com.typesafe.config.Config
import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine.db.ExecutionDatabaseKey
import cromwell.engine.io.IoInterface
import cromwell.engine.workflow.{CallKey, WorkflowOptions}
import cromwell.engine.{CallInputs, CallOutputs, HostInputs, _}
import cromwell.logging.WorkflowLogger
import cromwell.util.docker.SprayDockerRegistryApiClient
import org.slf4j.LoggerFactory
import wdl4s._

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

object Backend {
  class StdoutStderrException(message: String) extends RuntimeException(message)

  def from(backendConf: Config, actorSystem: ActorSystem): Backend = Backend.from(backendConf.getString("backend"), actorSystem)

  def from(name: String, actorSystem: ActorSystem) = {
    try {
      //      val constructor = Class.forName(name).getConstructor(ActorSystem.getClass)
      //      constructor.newInstance(actorSystem).asInstanceOf[Backend]
      Class.forName(name).newInstance().asInstanceOf[Backend]
    } catch {
      case exception: Exception =>
        throw new IllegalArgumentException(s"No backends could be loaded with the config $name. Reason: ${exception.getMessage}", exception)
    }
  }

  case class RestartableWorkflow(id: WorkflowId, source: WorkflowSourceFiles)
}

trait JobKey

/**
  * Trait to be implemented by concrete backends.
  */
trait Backend {
  type BackendCall <: backend.BackendCall

  def actorSystem: ActorSystem

  /**
    * Attempt to evaluate all the ${...} tags in a command and return a String representation
    * of the command.  This could fail for a variety of reasons related to expression evaluation
    * which is why it returns a Try[String]
    */
  def instantiateCommand(backendCall: BackendCall): Try[String] = {
    val backendInputs = adjustInputPaths(backendCall.key, backendCall.runtimeAttributes, backendCall.locallyQualifiedInputs, backendCall.workflowDescriptor)
    backendCall.call.instantiateCommandLine(backendInputs, backendCall.engineFunctions)
  }

  /**
    * Return a possibly altered copy of inputs reflecting any localization of input file paths that might have
    * been performed for this `Backend` implementation.
    */
  def adjustInputPaths(callKey: CallKey,
                       runtimeAttributes: RuntimeAttributes,
                       inputs: CallInputs,
                       workflowDescriptor: WorkflowDescriptor): CallInputs

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

  def workflowContext(workflowOptions: WorkflowOptions, workflowId: WorkflowId, name: String): WorkflowContext

  def engineFunctions(ioInterface: IoInterface, workflowContext: WorkflowContext): WorkflowEngineFunctions

  /**
    * Do any work that needs to be done <b>before</b> attempting to restart a workflow.
    */
  def prepareForRestart(restartableWorkflow: WorkflowDescriptor)(implicit ec: ExecutionContext): Future[Unit]

  /**
    * Return CallStandardOutput which contains the stdout/stderr of the particular call
    */
  def stdoutStderr(descriptor: WorkflowDescriptor, callName: String, index: ExecutionIndex): CallLogs

  def backendName: String

  /**
    * Validate that workflow options contain all required information.
    */
  @throws[IllegalArgumentException]("if a value is missing / incorrect")
  def assertWorkflowOptions(options: WorkflowOptions): Unit = {}

  private[backend] def backendClassString = backendName.toLowerCase.capitalize + "Backend"

  /** Default implementation assumes backends do not support resume, returns an empty Map. */
  def findResumableExecutions(id: WorkflowId)(implicit ec: ExecutionContext): Future[Map[ExecutionDatabaseKey, JobKey]] = Future.successful(Map.empty)

  def workflowLogger(descriptor: WorkflowDescriptor) = WorkflowLogger(
    backendClassString,
    descriptor,
    otherLoggers = Seq(LoggerFactory.getLogger(getClass.getName))
  )

  def workflowLoggerWithCall(backendCall: BackendCall) = WorkflowLogger(
    backendClassString,
    backendCall.workflowDescriptor,
    otherLoggers = Seq(LoggerFactory.getLogger(getClass.getName)),
    callTag = Option(backendCall.key.tag)
  )

  lazy val dockerHashClient = new SprayDockerRegistryApiClient()(actorSystem)
}
