package cromwell.engine.backend

import com.typesafe.config.{Config, ConfigFactory}
import cromwell.binding
import cromwell.binding.WdlExpression.ScopedLookupFunction
import cromwell.binding._
import cromwell.binding.types.WdlFileType
import cromwell.binding.values.{WdlFile, WdlString, WdlValue}
import cromwell.engine._
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.backend.local.LocalEngineFunctions
import cromwell.engine.backend.jes.JesBackend
import cromwell.engine.db.DataAccess
import cromwell.engine.backend.local.LocalBackend

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

object Backend {
  def from(backendConf: Config): Backend = {
    // FIXME: A log message would be a nice touch
    backendConf.getString("backend").toLowerCase match {
      case "local" => new LocalBackend
      case "jes" => new JesBackend
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

  /**
   * Presuming successful completion of the specified call, evaluate its outputs.
   */
  protected def evaluateCallOutputs(workflowDescriptor: WorkflowDescriptor,
                                    call: Call,
                                    hostAbsoluteFilePath: String => String,
                                    localEngineFunctions: LocalEngineFunctions,
                                    scopedLookupFunction: ScopedLookupFunction): Try[Map[String, WdlValue]] = {
    /**
     * Handle possible auto-conversion from an output expression `WdlString` to a `WdlFile` task output.
     * The following should work:
     *
     * <pre>
     * outputs {
     *   File bam = "foo.bam"
     * }
     * </pre>
     *
     * Output values that are not of type `WdlString` and are not being assigned to `WdlFiles` should be passed
     * through unaltered.  No other output conversions are currently supported and will result in `Failure`s.
     */
    def outputAutoConversion(callFqn: String, taskOutput: TaskOutput, rawOutputValue: WdlValue): Try[WdlValue] = {
      rawOutputValue match {
        // Autoconvert String -> File.
        case v: WdlString if taskOutput.wdlType == WdlFileType => Success(WdlFile(hostAbsoluteFilePath(v.value)))
        // Pass through matching types.
        case v if v.wdlType == taskOutput.wdlType => Success(v)
        // Fail other mismatched types.
        case _ =>
          val message = s"Expression '$rawOutputValue' of type ${rawOutputValue.wdlType.toWdlString} cannot be converted to ${taskOutput.wdlType.toWdlString} for output '$callFqn.${taskOutput.name}'."
          Failure(new RuntimeException(message))
      }
    }

    // Evaluate output expressions, performing conversions from String -> File where required.
    val outputMappings = call.task.outputs.map { taskOutput =>
      val tryConvertedValue =
        for {
          expressionValue <- taskOutput.expression.evaluate(scopedLookupFunction, localEngineFunctions)
          convertedValue <- outputAutoConversion(call.fullyQualifiedName, taskOutput, expressionValue)
        } yield convertedValue
      taskOutput.name -> tryConvertedValue
    }

    val taskOutputFailures = outputMappings.filter { _._2.isFailure }

    if (taskOutputFailures.isEmpty) {
      val unwrappedMap = outputMappings.collect { case (name, Success(wdlValue)) => name -> wdlValue }.toMap
      Success(unwrappedMap)
    } else {
      val message = taskOutputFailures.collect { case (name, Failure(e)) => s"$name: $e" }.mkString("\n")
      Failure(new Throwable(s"Workflow ${workflowDescriptor.id}: $message"))
    }
  }
}
