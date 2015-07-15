package cromwell.engine.backend.local

import java.io.{File, Writer}
import java.nio.file.{Files, Path, Paths}

import com.typesafe.scalalogging.LazyLogging
import cromwell.binding.WdlExpression.ScopedLookupFunction
import cromwell.binding._
import cromwell.binding.types.WdlFileType
import cromwell.binding.values.{WdlFile, WdlString, WdlValue}
import cromwell.engine._
import cromwell.engine.backend.Backend
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.db.{CallStatus, DataAccess}
import cromwell.util.FileUtil
import cromwell.util.FileUtil._

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.sys.process._
import scala.util.{Success, Failure, Try}


object LocalBackend {

  /**
   * Simple utility implicit class adding a method that writes a line and appends a newline in one shot.
   */
  implicit class WriteWithNewline(val writer: Writer) extends AnyVal {
    def writeWithNewline(string: String): Unit = {
      writer.write(string)
      writer.write("\n")
    }
  }
}


/**
 * Handles both local Docker runs as well as local direct command line executions.
 */
class LocalBackend extends Backend with LazyLogging {
  /**
   * Executes the specified command line, using the supplied lookup function for expression evaluation.
   * Returns a `Map[String, Try[WdlValue]]` of output names to values.
   */
  override def executeCommand(commandLine: String, 
                              workflowDescriptor: WorkflowDescriptor, 
                              call: Call, 
                              backendInputs: CallInputs,
                              scopedLookupFunction: ScopedLookupFunction): Try[Map[String, WdlValue]] = {
    import LocalBackend._

    val hostCallDirectory = Paths.get(hostExecutionPath(workflowDescriptor).toFile.getAbsolutePath, "call-" + call.name).toFile
    hostCallDirectory.mkdirs()

    val (stdoutFile, stdoutWriter) = FileUtil.tempFileAndWriter("stdout", hostCallDirectory)
    val (stderrFile, stderrWriter) = FileUtil.tempFileAndWriter("stderr", hostCallDirectory)
    val (commandFile, commandWriter) = FileUtil.tempFileAndWriter("command", hostCallDirectory)

    val hostWorkflowExecutionFile = hostExecutionPath(workflowDescriptor).toFile
    val containerExecutionDir = s"/root/${workflowDescriptor.id.toString}"

    val parentPath = if (call.docker.isDefined) Paths.get(containerExecutionDir).toString else hostWorkflowExecutionFile.getAbsolutePath
    val callDirectory = Paths.get(parentPath, s"call-${call.name}")
    commandWriter.writeWithNewline(s"cd $callDirectory")

    commandWriter.writeWithNewline(commandLine)
    commandWriter.flushAndClose()

    def buildDockerRunCommand(image: String): String =
      // -v maps the host workflow executions directory to /root/<workflow id> on the container.
      // -i makes the run interactive, required for the cat and <&0 shenanigans that follow.
      s"docker run -v ${hostWorkflowExecutionFile.getAbsolutePath}:$containerExecutionDir -i $image"

    // Build the docker run command if docker is defined in the RuntimeAttributes, otherwise just the empty string.
    val dockerRun = call.docker.map { buildDockerRunCommand }.getOrElse("")

    // The aforementioned shenanigans generate standard output and then a bash invocation that takes
    // commands from standard input.
    val argv = Seq("/bin/bash", "-c", s"cat $commandFile | $dockerRun /bin/bash <&0")
    logger.info("Executing call with argv: " + argv)

    // The ! to the ProcessLogger captures standard output and error.
    val rc: Int = argv ! ProcessLogger(stdoutWriter writeWithNewline, stderrWriter writeWithNewline)
    Vector(stdoutWriter, stderrWriter).foreach {
      _.flushAndClose()
    }

    /**
     * Return a host absolute file path.
     */
    def hostAbsoluteFilePath(pathString: String): String =
      if (new File(pathString).isAbsolute) pathString else Paths.get(hostCallDirectory.toString, pathString).toString

    /**
     * Handle possible auto-conversion from string literal to WdlFile.
     */
    def possiblyAutoConvertedValue(taskOutput: TaskOutput, wdlValue: WdlValue): WdlValue = {
      wdlValue match {
        case v: WdlString =>
          taskOutput.wdlType match {
            case WdlFileType =>
              WdlFile(hostAbsoluteFilePath(v.value))
            case _ => v
          }
        case v => v
      }
    }

    val outputMappings = call.task.outputs.map { taskOutput =>
      val rawValue = taskOutput.expression.evaluate(
        scopedLookupFunction,
        new LocalEngineFunctions(TaskExecutionContext(stdoutFile, stderrFile))
      )
      taskOutput.name -> rawValue.map {
        possiblyAutoConvertedValue(taskOutput, _)
      }
    }

    if (rc == 0) {
      val taskOutputEvaluationFailures = outputMappings.filter { _._2.isFailure }

      if (taskOutputEvaluationFailures.isEmpty) {
        val unwrappedMap = outputMappings.collect { case (name, Success(wdlValue) ) => name -> wdlValue }.toMap
        Success(unwrappedMap)
      } else {
        val message = taskOutputEvaluationFailures.collect { case (name, Failure(e))  => s"$name: $e" }.mkString("\n")
        Failure(new Throwable(s"Workflow ${workflowDescriptor.id}: $message"))
      }
    } else {
      Failure(new Throwable(s"Workflow ${workflowDescriptor.id}: return code $rc for command: $commandLine"))
    }
  }

  /**
   * Creates host execution directory, inputs path, and outputs path.  Stages any input files into the workflow-inputs
   * directory and localizes their paths relative to the container.
   */
  override def initializeForWorkflow(descriptor: WorkflowDescriptor): HostInputs = {
    val hostExecutionDirectory = hostExecutionPath(descriptor).toFile
    hostExecutionDirectory.mkdirs()
    // FIXME: I don't think we need to do the next line, can't we do Paths.resolve on it 2 lines down?
    val hostExecutionAbsolutePath = hostExecutionDirectory.getAbsolutePath
    Array("workflow-inputs", "workflow-outputs") foreach { Paths.get(hostExecutionAbsolutePath, _).toFile.mkdir() }
    stageWorkflowInputs(descriptor)
  }

  /**
   * {{{cromwell-executions + workflow.name + workflow.id = cromwell-executions/three-step/0f00-ba4}}}
   */
  private def hostExecutionPath(workflow: WorkflowDescriptor): Path = Paths.get("cromwell-executions", workflow.name, workflow.id.toString)

  /**
   * Given the specified workflow descriptor and inputs, stage any WdlFiles into the workflow-inputs subdirectory
   * of the workflow execution directory.  Return a Map of the input values with any input WdlFiles adjusted to
   * reflect host paths.
   *
   * Original input path: /could/be/anywhere/input.bam
   * Host inputs path: $PWD/cromwell-executions/some-workflow-name/0f00-ba4/workflow-inputs/input.bam
   */
  private def stageWorkflowInputs(descriptor: WorkflowDescriptor): HostInputs = {

    val hostInputsPath = Paths.get(hostExecutionPath(descriptor).toFile.getAbsolutePath, "workflow-inputs")

    def stageInput(nameAndValue: (String, WdlValue)): (String, WdlValue) = {
      val (name, value) = nameAndValue
      val hostPathAdjustedValue = value match {
        case WdlFile(originalPath) =>
          val executionPath = Paths.get(hostInputsPath.toFile.getAbsolutePath, Paths.get(originalPath).getFileName.toString)
          Files.copy(Paths.get(originalPath), executionPath)
          WdlFile(executionPath.toString)
        case x => x
      }
      name -> hostPathAdjustedValue
    }
    descriptor.actualInputs map stageInput
  }

  /**
   * Return a possibly altered copy of inputs reflecting any localization of input file paths that might have
   * been performed for this `Backend` implementation.
   */
  override def adjustInputPaths(call: Call, inputs: CallInputs): CallInputs = {
    // If this call is using Docker, adjust input paths, otherwise return unaltered input paths.
    val CromwellExecutions = "cromwell-executions"
    def adjustPath(nameAndValue: (String, WdlValue)): (String, WdlValue) = {
      val (name, value) = nameAndValue
      val adjusted = value match {
        case WdlFile(path) =>
          // Host path would look like cromwell-executions/three-step/f00ba4/call-ps/stdout.txt
          // Container path should look like /root/f00ba4/call-ps/stdout.txt
          val fullPath = Paths.get(path).toFile.getAbsolutePath
          // Strip out everything before cromwell-executions.
          val pathUnderCromwellExecutions = fullPath.substring(fullPath.indexOf(CromwellExecutions) + CromwellExecutions.length)
          // Strip out the workflow name (the first component under cromwell-executions).
          val pathWithWorkflowName = Paths.get(pathUnderCromwellExecutions)
          WdlFile(WdlFile.appendPathsWithSlashSeparators("/root", pathWithWorkflowName.subpath(1, pathWithWorkflowName.getNameCount).toString))
        case x => x
      }
      name -> adjusted
    }

    // If this call is using Docker, adjust input paths, otherwise return unaltered input paths.
    if (call.docker.isDefined) inputs map adjustPath else inputs
  }

  override def adjustOutputPaths(call: Call, outputs: CallOutputs): CallOutputs = outputs

  /**
   * LocalBackend needs to force non-terminal calls back to NotStarted on restart.
   */
  override def handleCallRestarts(restartableWorkflows: Seq[RestartableWorkflow], dataAccess: DataAccess)
                                 (implicit ec: ExecutionContext): Future[Any] = {
    // Remove terminal states and the NotStarted state from the states which need to be reset to NotStarted.
    val StatusesNeedingUpdate = ExecutionStatus.values -- Set(ExecutionStatus.Failed, ExecutionStatus.Done, ExecutionStatus.NotStarted)
    def updateNonTerminalCalls(workflowId: WorkflowId, callFqnsToStatuses: Map[FullyQualifiedName, CallStatus]): Future[Unit] = {
      val callFqnsNeedingUpdate = callFqnsToStatuses.collect { case (callFqn, status) if StatusesNeedingUpdate.contains(status) => callFqn}
      dataAccess.setStatus(workflowId, callFqnsNeedingUpdate, ExecutionStatus.NotStarted)
    }

    val seqOfFutures = restartableWorkflows map { workflow =>
      for {
        callsToStatuses <- dataAccess.getExecutionStatuses(workflow.id)
        _ <- updateNonTerminalCalls(workflow.id, callsToStatuses)
      } yield ()
    }
    // The caller doesn't care about the result value, only that it's a Future.
    Future.sequence(seqOfFutures)
  }
}
