package cromwell.engine.backend.local

import java.io.{File, Writer}
import java.nio.file.{Files, Path, Paths}
import java.util.UUID

import com.typesafe.scalalogging.LazyLogging
import cromwell.binding.WdlExpression.ScopedLookupFunction
import cromwell.binding._
import cromwell.binding.values.{WdlArray, WdlFile, WdlValue}
import cromwell.engine.backend.{StdoutStderr, Backend}
import cromwell.engine.backend.Backend.{RestartableWorkflow, StdoutStderrException}
import cromwell.engine.db.{CallStatus, DataAccess}
import cromwell.engine.{ExecutionStatus, WorkflowId}
import cromwell.parser.BackendType
import cromwell.util.FileUtil
import cromwell.util.FileUtil._
import org.apache.commons.io.FileUtils

import scala.collection.JavaConversions._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.sys.process._
import scala.util.{Failure, Success, Try}

object LocalBackend {

  val CromwellExecutions = "cromwell-executions"

  /**
   * Simple utility implicit class adding a method that writes a line and appends a newline in one shot.
   */
  implicit class WriteWithNewline(val writer: Writer) extends AnyVal {
    def writeWithNewline(string: String): Unit = {
      writer.write(string)
      writer.write("\n")
    }
  }

  def findTempFile(root: Path, prefix: String) = {
    val regex = s"$prefix.*\\.tmp$$".r.unanchored
    val filesWithPrefix = Try(Files.newDirectoryStream(root)).map(
      stream => stream.iterator().toIterator.toList.collect {
        case path:Path if regex.findFirstIn(path.toAbsolutePath.toString).isDefined => path
      }
    )
    filesWithPrefix match {
      case Success(paths) => paths match {
        case head :: Nil => WdlFile(head.toAbsolutePath.toString)
        case Nil => throw new StdoutStderrException(s"No $prefix file found")
        case s => throw new StdoutStderrException(s"Multiple files matched with prefix $prefix:\n${s.map(_.toString).mkString(", ")}")
      }
      case Failure(ex) => throw ex
    }
  }

  /**
   * {{{cromwell-executions + workflow.name + workflow.id = cromwell-executions/three-step/0f00-ba4}}}
   */
  def hostExecutionPath(workflow: WorkflowDescriptor): Path =
    hostExecutionPath(workflow.name, workflow.id)

  def hostExecutionPath(workflowName: String, workflowUuid: UUID): Path =
    Paths.get(CromwellExecutions, workflowName, workflowUuid.toString)

  def hostCallPath(workflow: WorkflowDescriptor, callName: String): Path =
    Paths.get(hostExecutionPath(workflow).toFile.getAbsolutePath, s"call-$callName")

  def hostCallPath(workflowName: String, workflowUuid: UUID, callName: String): Path =
    Paths.get(hostExecutionPath(workflowName, workflowUuid).toFile.getAbsolutePath, s"call-$callName")
}


/**
 * Handles both local Docker runs as well as local direct command line executions.
 */
class LocalBackend extends Backend with LazyLogging {

  import LocalBackend._

  override def stdoutStderr(workflowId: WorkflowId, workflowName: String, callName: String): StdoutStderr = {
    val dir = hostCallPath(workflowName, workflowId, callName)
    StdoutStderr(
      stdout = findTempFile(dir, prefix = "stdout"),
      stderr = findTempFile(dir, prefix = "stderr")
    )
  }

  /**
   * Executes the specified command line, using the supplied lookup function for expression evaluation.
   * Returns a `Map[String, Try[WdlValue]]` of output names to values.
   */
  override def executeCommand(instantiatedCommandLine: String, 
                              workflowDescriptor: WorkflowDescriptor, 
                              call: Call, 
                              backendInputs: CallInputs,
                              scopedLookupFunction: ScopedLookupFunction): Try[Map[String, WdlValue]] = {

    val workflowRootAbsolutePathOnHost = hostExecutionPath(workflowDescriptor).toFile.getAbsolutePath

    val hostCallDirectory = hostCallPath(workflowDescriptor, call.name).toFile
    hostCallDirectory.mkdirs()

    val (stdoutFile, stdoutWriter) = FileUtil.tempFileAndWriter("stdout", hostCallDirectory)
    val (stderrFile, stderrWriter) = FileUtil.tempFileAndWriter("stderr", hostCallDirectory)
    val (commandFile, commandWriter) = FileUtil.tempFileAndWriter("command", hostCallDirectory)

    val dockerContainerExecutionDir = s"/root/${workflowDescriptor.id.toString}"

    val parentPath = if (call.docker.isDefined) Paths.get(dockerContainerExecutionDir).toString else workflowRootAbsolutePathOnHost
    val callDirectory = Paths.get(parentPath, s"call-${call.name}")

    commandWriter.writeWithNewline(s"cd $callDirectory")
    commandWriter.writeWithNewline(instantiatedCommandLine)
    commandWriter.flushAndClose()

    def buildDockerRunCommand(image: String): String =
      // -v maps the host workflow executions directory to /root/<workflow id> on the container.
      // -i makes the run interactive, required for the cat and <&0 shenanigans that follow.
      s"docker run -v $workflowRootAbsolutePathOnHost:$dockerContainerExecutionDir -i $image"

    // Build the docker run command if docker is defined in the RuntimeAttributes, otherwise just the empty string.
    val dockerRun = call.docker.map { buildDockerRunCommand }.getOrElse("")

    // The aforementioned shenanigans generate standard output and then a bash invocation that takes
    // commands from standard input.
    val argv = Seq("/bin/bash", "-c", s"cat $commandFile | $dockerRun /bin/bash <&0")
    logger.debug("Executing call with argv: " + argv)

    // The ! to the ProcessLogger captures standard output and error.
    val rc: Int = argv ! ProcessLogger(stdoutWriter writeWithNewline, stderrWriter writeWithNewline)
    Vector(stdoutWriter, stderrWriter).foreach { _.flushAndClose() }

    /**
     * Return a host absolute file path.
     */
    def hostAbsoluteFilePath(pathString: String): String =
      if (new File(pathString).isAbsolute) pathString else Paths.get(hostCallDirectory.toString, pathString).toString

    val localEngineFunctions = new LocalEngineFunctions(TaskExecutionContext(stdoutFile, stderrFile, Paths.get(hostCallDirectory.getAbsolutePath)))

    val stderrFileLength = new File(stderrFile.toString).length

    if (call.failOnStderr && stderrFileLength > 0) {
      Failure(new Throwable(s"Workflow ${workflowDescriptor.id}: stderr has length $stderrFileLength for command: $instantiatedCommandLine"))
    } else {
      if (rc == 0) {
        evaluateCallOutputs(workflowDescriptor, call, hostAbsoluteFilePath, localEngineFunctions, scopedLookupFunction, interpolateStrings=true)
      } else {
        Failure(new Throwable(s"Workflow ${workflowDescriptor.id}: return code $rc for command: $instantiatedCommandLine\n\nFull command was: ${argv.mkString(" ")}"))
      }
    }
  }

  /**
   * Creates host execution directory, inputs path, and outputs path.  Stages any input files into the workflow-inputs
   * directory and localizes their paths relative to the container.
   */
  override def initializeForWorkflow(descriptor: WorkflowDescriptor): HostInputs = {
    val hostExecutionDirectory = hostExecutionPath(descriptor).toFile
    hostExecutionDirectory.mkdirs()
    val hostExecutionAbsolutePath = hostExecutionDirectory.getAbsolutePath
    Array("workflow-inputs", "workflow-outputs") foreach { Paths.get(hostExecutionAbsolutePath, _).toFile.mkdir() }
    stageWorkflowInputs(descriptor)
  }

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
    descriptor.actualInputs map {case(name, value) => name -> stageWdlValue(value, hostInputsPath)}
  }

  private def stageInputArray(array: WdlArray, hostInputsPath: Path): WdlArray = array.map(stageWdlValue(_, hostInputsPath))

  private def stageWdlValue(value: WdlValue, hostInputsPath: Path): WdlValue = value match {
    case w:WdlFile => stageWdlFile(w, hostInputsPath)
    case a:WdlArray => stageInputArray(a, hostInputsPath)
    case x => x
  }

  private def stageWdlFile(wdlFile: WdlFile, hostInputsPath: Path): WdlFile = {
    val originalPath = Paths.get(wdlFile.value)
    val executionPath = hostInputsPath.resolve(originalPath.getFileName.toString)
    if (Files.isDirectory(originalPath)) {
      FileUtils.copyDirectory(originalPath.toFile, executionPath.toFile)
    } else {
      FileUtils.copyFile(originalPath.toFile, executionPath.toFile)
    }
    WdlFile(executionPath.toString)
  }

  /**
   * Return a possibly altered copy of inputs reflecting any localization of input file paths that might have
   * been performed for this `Backend` implementation.
   */
  override def adjustInputPaths(call: Call, inputs: CallInputs): CallInputs = {
    // If this call is using Docker, adjust input paths, otherwise return unaltered input paths.
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

  override def backendType = BackendType.LOCAL
}
