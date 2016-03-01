package cromwell.engine.backend

import java.nio.file.Path

import akka.actor.ActorSystem
import com.google.api.client.util.ExponentialBackOff
import com.typesafe.config.Config
import cromwell.engine._
import cromwell.engine.backend.jes.JesBackend
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.backend.runtimeattributes.CromwellRuntimeAttributes
import cromwell.engine.backend.sge.SgeBackend
import cromwell.engine.db.DataAccess.ExecutionKeyToJobKey
import cromwell.engine.db.ExecutionDatabaseKey
import cromwell.engine.io.IoInterface
import cromwell.engine.workflow.{CallKey, WorkflowOptions}
import cromwell.engine.workflow.{BackendCallKey, WorkflowOptions}
import cromwell.engine.{HostInputs, CallOutputs}
import cromwell.logging.WorkflowLogger
import cromwell.util.docker.SprayDockerRegistryApiClient
import org.slf4j.LoggerFactory
import wdl4s._
import wdl4s.values.WdlValue

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Success, Try}

object Backend {
  class StdoutStderrException(message: String) extends RuntimeException(message)

  def from(backendConf: Config, actorSystem: ActorSystem): Backend = Backend.from(backendConf.getString("backend"), actorSystem)

  def from(backendType: BackendType, actorSystem: ActorSystem): Backend = backendType match {
    case BackendType.LOCAL => LocalBackend(actorSystem)
    case BackendType.JES => JesBackend(actorSystem)
    case BackendType.SGE => SgeBackend(actorSystem)
    case doh => throw new IllegalArgumentException(s"$doh is not a recognized backend")
  }

  def from(name: String, actorSystem: ActorSystem): Backend = {
    val backendType = name.toBackendType
    Backend.from(backendType, actorSystem)
  }

  case class RestartableWorkflow(id: WorkflowId, source: WorkflowSourceFiles)

  implicit class BackendyString(val backendType: String) extends AnyVal {
    def toBackendType: BackendType = {
      try {
        BackendType.valueOf(backendType.toUpperCase)
      } catch {
        case e: Exception => throw new IllegalArgumentException(s"$backendType is not a recognized backend")
      }
    }
  }
}

trait JobKey

final case class AttemptedLookupResult(name: String, value: Try[WdlValue]) {
  def toPair = name -> value
}

object AttemptedLookupResult {
  implicit class AugmentedAttemptedLookupSequence(s: Seq[AttemptedLookupResult]) {
    def toLookupMap: Map[String, WdlValue] = s collect {
      case AttemptedLookupResult(name, Success(value)) => (name, value)
    } toMap
  }
}

/**
 * Trait to be implemented by concrete backends.
 */
trait Backend {
  type BackendCall <: backend.BackendCall

  def actorSystem: ActorSystem

  def rootPath(workflowOptions: WorkflowOptions): String

  /**
   * Return a possibly altered copy of inputs reflecting any localization of input file paths that might have
   * been performed for this `Backend` implementation.
   */
  def adjustInputPaths(backendCall: BackendCall): CallInputs

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
               key: BackendCallKey,
               locallyQualifiedInputs: CallInputs = Map.empty[String, WdlValue],
               abortRegistrationFunction: Option[AbortRegistrationFunction] = None): BackendCall

  def engineFunctions(ioInterface: IoInterface, workflowContext: WorkflowContext): WorkflowEngineFunctions

  /**
   * Do any work that needs to be done <b>before</b> attempting to restart a workflow.
   */
  def prepareForRestart(restartableWorkflow: WorkflowDescriptor)(implicit ec: ExecutionContext): Future[Unit]

  /**
   * Return CallLogs which contains the stdout/stderr of the particular call
   */
  def stdoutStderr(backendCall: BackendCall): CallLogs

  def backendType: BackendType

  /**
   * Validate that workflow options contain all required information.
   */
  @throws[IllegalArgumentException]("if a value is missing / incorrect")
  def assertWorkflowOptions(options: WorkflowOptions): Unit = {}

  private[backend] def backendClassString = backendType.toString.toLowerCase.capitalize + "Backend"

  /** Default implementation assumes backends do not support resume, returns an empty Map. */
  def findResumableExecutions(id: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[ExecutionKeyToJobKey]] = Future.successful(List.empty)

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

  def pollBackoff: ExponentialBackOff

  def executionInfoKeys: List[String]
}
