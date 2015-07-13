package cromwell.engine.backend

import com.typesafe.config.{Config, ConfigFactory}
import cromwell.binding
import cromwell.binding.WdlExpression.ScopedLookupFunction
import cromwell.binding._
import cromwell.binding.values.WdlValue
import cromwell.engine._
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.backend.jes.JesBackend
import cromwell.engine.db.DataAccess
import cromwell.engine.backend.local.LocalBackend

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

object Backend {
  def from(backendConf: Config): Backend = {
    backendConf.getString("backend").toLowerCase match {
      case "local" =>
        println("FOO: LOCAL BACKEND")
        new LocalBackend
      case "jes" =>
        println("FOO: JES BACKEND")
        throw new Exception("BARF")
        new JesBackend
      case doh => throw new IllegalArgumentException(s"$doh is not a recognized backend")
    }
  }

  case class RestartableWorkflow(id: WorkflowId, source: WdlSource, json: WdlJson, inputs: binding.WorkflowRawInputs)
}

/**
 * Trait to be implemented by concrete backends.
 */
trait Backend {

  /**
   * Return a possibly altered copy of inputs reflecting any localization of input file paths that might have
   * been performed for this `Backend` implementation.
   */
  def adjustInputPaths(call: Call, inputs: CallInputs): CallInputs

  /**
   * Do whatever work is required to initialize the workflow, returning a copy of
   * the coerced inputs present in the `WorkflowDescriptor` with any input `WdlFile`s
   * adjusted for the host workflow execution path.
   */
  def initializeForWorkflow(workflow: WorkflowDescriptor): HostInputs

  /**
   * Execute the specified command line using the provided symbol store, evaluating the task outputs to produce
   * a mapping of local task output names to WDL values.
   */
  def executeCommand(commandLine: String, 
                     workflowDescriptor: WorkflowDescriptor, 
                     call: Call, 
                     backendInputs: CallInputs, 
                     scopedLookupFunction: ScopedLookupFunction): Try[Map[String, WdlValue]]

  /**
   * Do whatever is appropriate for this backend implementation to support restarting the specified workflows.
   */
  def handleCallRestarts(restartableWorkflows: Seq[RestartableWorkflow], dataAccess: DataAccess)(implicit ec: ExecutionContext): Future[Any]
}
