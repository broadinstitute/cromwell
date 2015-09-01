package cromwell.engine.backend.local

import java.io.Writer
import java.nio.file.{Files, Path, Paths}

import com.typesafe.scalalogging.LazyLogging
import cromwell.binding._
import cromwell.engine._
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.backend.{Backend, TaskAbortedException}
import cromwell.engine.db.{CallStatus, DataAccess}
import cromwell.parser.BackendType
import cromwell.util.FileUtil._

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
  type BackendCall = LocalBackendCall

  import LocalBackend._

  override def bindCall(workflowDescriptor: WorkflowDescriptor,
                        call: Call,
                        locallyQualifiedInputs: CallInputs,
                        abortRegistrationFunction: AbortRegistrationFunction): BackendCall = {
    LocalBackendCall(this, workflowDescriptor, call, locallyQualifiedInputs, abortRegistrationFunction)
  }

  override def execute(backendCall: BackendCall): Try[CallOutputs] =  {
    val tag = makeTag(backendCall)
    backendCall.instantiateCommand match {
      case Success(instantiatedCommand) =>
        logger.info(s"$tag `$instantiatedCommand`")
        writeScript(backendCall, instantiatedCommand, backendCall.containerCallRoot)
        runSubprocess(backendCall)
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

  /**
   * Writes the script file containing the user's command from the WDL as well
   * as some extra shell code for monitoring jobs
   */
  private def writeScript(backendCall: BackendCall, instantiatedCommand: String, containerRoot: Path) = {
    val (_, scriptWriter) = backendCall.script.fileAndWriter
    scriptWriter.writeWithNewline("#!/bin/sh")
    scriptWriter.writeWithNewline(s"cd $containerRoot")
    scriptWriter.writeWithNewline(instantiatedCommand)
    scriptWriter.writeWithNewline("echo $? > rc")
    scriptWriter.flushAndClose()
  }

  /**
   * --rm automatically deletes the container upon exit
   * -v maps the host workflow executions directory to /root/<workflow id> on the container.
   * -i makes the run interactive, required for the cat and <&0 shenanigans that follow.
   */
  private def buildDockerRunCommand(backendCall: BackendCall, image: String): String =
    s"docker run --rm -v ${backendCall.workflowRootPath.toAbsolutePath}:${backendCall.dockerContainerExecutionDir} -i $image"

  private def runSubprocess(backendCall: BackendCall): Try[CallOutputs] = {
    val tag = makeTag(backendCall)
    val (_, stdoutWriter) = backendCall.stdout.fileAndWriter
    val (_, stderrWriter) = backendCall.stderr.fileAndWriter
    val dockerRun = backendCall.call.docker.map(d => buildDockerRunCommand(backendCall, d)).getOrElse("")
    val argv = Seq("/bin/bash", "-c", s"cat ${backendCall.script} | $dockerRun /bin/bash <&0")
    val process = argv.run(ProcessLogger(stdoutWriter writeWithNewline, stderrWriter writeWithNewline))

    // TODO: As currently implemented, this process.destroy() will kill the bash process but *not* its descendants. See ticket DSDEEPB-848.
    backendCall.callAbortRegistrationFunction.register(AbortFunction(() => process.destroy()))
    val backendCommandString = argv.map(s => "\""+s+"\"").mkString(" ")
    logger.info(s"$tag command: $backendCommandString")
    process.exitValue() // blocks until process finishes
    Vector(stdoutWriter, stderrWriter) foreach { _.flushAndClose() }

    val stderrFileLength = Files.size(backendCall.stderr)
    val rc = Try(backendCall.rc.slurp.stripLineEnd.toInt)
    if (backendCall.call.failOnStderr && stderrFileLength > 0) {
      Failure(new Throwable(s"Call ${backendCall.call.fullyQualifiedName}, Workflow ${backendCall.workflowDescriptor.id}: stderr has length $stderrFileLength"))
    } else {
      rc match {
        case Success(0) => postProcess(backendCall)
        case Success(143) => Failure(new TaskAbortedException()) // Special case to check for SIGTERM exit code - implying abort
        case _ => Failure(new Throwable(s"Call ${backendCall.call.fullyQualifiedName}, Workflow ${backendCall.workflowDescriptor.id}: $rc\n\nFull command was: ${argv.mkString(" ")}"))
      }
    }
  }
}
