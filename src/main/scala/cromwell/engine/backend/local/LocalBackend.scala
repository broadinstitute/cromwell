package cromwell.engine.backend.local

import java.io.{File, Writer}
import java.nio.file.{Files, Path, Paths}

import com.typesafe.scalalogging.LazyLogging
import cromwell.binding._
import cromwell.binding.types.{WdlFileType, WdlMapType}
import cromwell.binding.values._
import cromwell.engine._
import cromwell.engine.backend.Backend.{RestartableWorkflow, StdoutStderrException}
import cromwell.engine.backend.{Backend, BackendCall, StdoutStderr, TaskAbortedException}
import cromwell.engine.db.{CallStatus, DataAccess}
import cromwell.parser.BackendType
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

  /**
   * {{{cromwell-executions + workflow.name + workflow.id = cromwell-executions/three-step/0f00-ba4}}}
   */
  def hostExecutionPath(workflow: WorkflowDescriptor): Path =
    hostExecutionPath(workflow.name, workflow.id)

  def hostExecutionPath(workflowName: String, workflowUuid: WorkflowId): Path =
    Paths.get(CromwellExecutions, workflowName, workflowUuid.id.toString)

  def hostCallPath(workflow: WorkflowDescriptor, callName: String): Path =
    Paths.get(hostExecutionPath(workflow).toFile.getAbsolutePath, s"call-$callName")

  def hostCallPath(workflowName: String, workflowUuid: WorkflowId, callName: String): Path =
    Paths.get(hostExecutionPath(workflowName, workflowUuid).toFile.getAbsolutePath, s"call-$callName")
}


/**
 * Handles both local Docker runs as well as local direct command line executions.
 */
class LocalBackend extends Backend with LocalFileSystemOperations with LazyLogging {

  import LocalBackend._

  override def bindCall(workflowDescriptor: WorkflowDescriptor,
                        call: Call,
                        locallyQualifiedInputs: CallInputs,
                        abortRegistrationFunction: AbortFunctionRegistration): BackendCall = {
    LocalBackendCall(this, workflowDescriptor, call, locallyQualifiedInputs, abortRegistrationFunction)
  }

  override def execute(backendCall: BackendCall): Try[CallOutputs] = backendCall match {
    case bc: LocalBackendCall => execute(bc)
    case _ => Failure(new UnsupportedOperationException("LocalBackend can only execute LocalBackendCalls"))
  }

  /**
   * Executes the specified command line, using the supplied lookup function for expression evaluation.
   * Returns a `Map[CallOutputs]` of output names to values.
   */
  private def execute(backendCall: LocalBackendCall): Try[CallOutputs]  = {
    val tag: String = s"${this.getClass.getName} [UUID(${backendCall.workflowDescriptor.shortId}):${backendCall.call.name}]"
    val dockerContainerExecutionDir = s"/root/${backendCall.workflowDescriptor.id.toString}"
    val containerCallRoot = backendCall.call.docker match {
      case Some(docker) => Paths.get(dockerContainerExecutionDir).resolve(s"call-${backendCall.call.name}")
      case None => backendCall.callRootPath
    }
    def buildDockerRunCommand(image: String): String =
      // -v maps the host workflow executions directory to /root/<workflow id> on the container.
      // -i makes the run interactive, required for the cat and <&0 shenanigans that follow.
      s"docker run -v ${backendCall.workflowRootPath.toAbsolutePath}:$dockerContainerExecutionDir -i $image"

    backendCall.instantiateCommand match {
      case Success(instantiatedCommand) =>
        logger.info(s"$tag `$instantiatedCommand`")

        val (stdoutFile, stdoutWriter) = backendCall.stdout.fileAndWriter
        val (stderrFile, stderrWriter) = backendCall.stderr.fileAndWriter
        val (scriptFile, scriptWriter) = backendCall.script.fileAndWriter
        scriptWriter.writeWithNewline(s"cd $containerCallRoot")
        scriptWriter.writeWithNewline(s"$instantiatedCommand")
        scriptWriter.flushAndClose()

        val dockerRun = backendCall.call.docker.map(buildDockerRunCommand).getOrElse("")
        val argv = Seq("/bin/bash", "-c", s"cat $scriptFile | $dockerRun /bin/bash <&0")
        val backendCommandString = argv.map(s => "\""+s+"\"").mkString(" ")
        logger.debug(s"$tag backend command: $backendCommandString")
        val process = argv.run(ProcessLogger(stdoutWriter writeWithNewline, stderrWriter writeWithNewline))

        // TODO: As currently implemented, this process.destroy() will kill the bash process but *not* its descendants. See ticket DSDEEPB-848.
        backendCall.callAbortRegistrationFunction.register(AbortFunction(() => process.destroy()))
        val rc: Int = process.exitValue()
        Vector(stdoutWriter, stderrWriter) foreach { _.flushAndClose() }

        val stderrFileLength = new File(stderrFile.toString).length

        if (backendCall.call.failOnStderr && stderrFileLength > 0) {
          Failure(new Throwable(s"Workflow ${backendCall.workflowDescriptor.id}: stderr has length $stderrFileLength for command: $instantiatedCommand"))
        } else {
          rc match {
            case 0 => postProcess(backendCall)
            case 143 => Failure(new TaskAbortedException()) // Special case to check for SIGTERM exit code - implying abort
            case _ => Failure(new Throwable(s"Workflow ${backendCall.workflowDescriptor.id}: return code $rc for command: $instantiatedCommand\n\nFull command was: ${argv.mkString(" ")}"))
          }
        }
      case Failure(ex) => Failure(ex)
    }
  }


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
