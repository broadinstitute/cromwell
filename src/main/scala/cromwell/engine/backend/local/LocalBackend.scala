package cromwell.engine.backend.local

import java.io.Writer
import java.nio.file.{Files, Path, Paths}

import akka.actor.ActorSystem
import better.files._
import cromwell.binding._
import cromwell.engine.ExecutionIndex._
import cromwell.engine._
import cromwell.engine.backend._
import cromwell.engine.db.DataAccess._
import cromwell.engine.db.{CallStatus, ExecutionDatabaseKey}
import cromwell.engine.workflow.CallKey
import cromwell.parser.BackendType
import cromwell.util.FileUtil._
import org.joda.time.DateTime

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.sys.process._
import scala.util.{Failure, Success, Try}

object LocalBackend {

  val ContainerRoot = "/root"
  val CallPrefix = "call"
  val ShardPrefix = "shard"

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
    Paths.get(SharedFileSystem.CromwellExecutionRoot, workflowName, workflowUuid.id.toString)

  def hostCallPath(workflow: WorkflowDescriptor, callName: String, callIndex: ExecutionIndex): Path = {
    hostCallPath(workflow.name, workflow.id, callName, callIndex)
  }

  def hostCallPath(workflowName: String, workflowUuid: WorkflowId, callName: String, callIndex: ExecutionIndex): Path =  {
    val rootCallPath = hostExecutionPath(workflowName, workflowUuid).resolve(s"$CallPrefix-$callName")
    callIndex match {
      case Some(index) => rootCallPath.resolve(s"$ShardPrefix-$index")
      case None => rootCallPath
    }
  }

  /**
   * Root workflow execution path for container.
   */
  def containerExecutionPath(workflow: WorkflowDescriptor): Path = Paths.get(ContainerRoot, workflow.id.toString)

  /**
   * Root Call execution path for container.
   */
  def containerCallPath(workflow: WorkflowDescriptor, callName: String, callIndex: ExecutionIndex): Path = {
    val rootCallPath = containerExecutionPath(workflow).resolve(s"$CallPrefix-$callName")
    callIndex match {
      case Some(index) => rootCallPath.resolve(s"$ShardPrefix-$index")
      case None => rootCallPath
    }
  }
}


/**
 * Handles both local Docker runs as well as local direct command line executions.
 */
case class LocalBackend(actorSystem: ActorSystem) extends Backend with SharedFileSystem {
  type BackendCall = LocalBackendCall

  import LocalBackend._

  override def bindCall(workflowDescriptor: WorkflowDescriptor,
                        key: CallKey,
                        locallyQualifiedInputs: CallInputs,
                        abortRegistrationFunction: AbortRegistrationFunction): BackendCall = {
    LocalBackendCall(this, workflowDescriptor, key, locallyQualifiedInputs, abortRegistrationFunction)
  }

  def execute(backendCall: BackendCall)(implicit ec: ExecutionContext): Future[ExecutionHandle] = Future({
    val logger = workflowLoggerWithCall(backendCall)
    backendCall.instantiateCommand match {
      case Success(instantiatedCommand) =>
        logger.info(s"`$instantiatedCommand`")
        writeScript(backendCall, instantiatedCommand, backendCall.containerCallRoot)
        runSubprocess(backendCall)
      case Failure(ex) => FailedExecution(ex).future
    }
  }).flatten map CompletedExecutionHandle

  /**
   * LocalBackend needs to force non-terminal calls back to NotStarted on restart.
   */
  override def prepareForRestart(restartableWorkflow: WorkflowDescriptor)(implicit ec: ExecutionContext): Future[Unit] = {
    // Remove terminal states and the NotStarted state from the states which need to be reset to NotStarted.
    val StatusesNeedingUpdate = ExecutionStatus.values -- Set(ExecutionStatus.Failed, ExecutionStatus.Done, ExecutionStatus.NotStarted)
    def updateNonTerminalCalls(workflowId: WorkflowId, keyToStatusMap: Map[ExecutionDatabaseKey, CallStatus]): Future[Unit] = {
      val callFqnsNeedingUpdate = keyToStatusMap collect { case (callFqn, callStatus) if StatusesNeedingUpdate.contains(callStatus.executionStatus) => callFqn }
      globalDataAccess.setStatus(workflowId, callFqnsNeedingUpdate, CallStatus(ExecutionStatus.NotStarted, None, None, None))
    }

    for {
      callsToStatuses <- globalDataAccess.getExecutionStatuses(restartableWorkflow.id)
      _ <- updateNonTerminalCalls(restartableWorkflow.id, callsToStatuses)
    } yield ()
  }

  override def backendType = BackendType.LOCAL

  /**
   * Writes the script file containing the user's command from the WDL as well
   * as some extra shell code for monitoring jobs
   */
  private def writeScript(backendCall: BackendCall, instantiatedCommand: String, containerRoot: Path) = {
    backendCall.script.write(
      s"""#!/bin/sh
         |cd $containerRoot
         |$instantiatedCommand
         |echo $$? > rc
         |""".stripMargin)
  }

  /**
   * --rm automatically deletes the container upon exit
   * -v maps the host workflow executions directory to /root/<workflow id> on the container.
   * -i makes the run interactive, required for the cat and <&0 shenanigans that follow.
   */
  private def buildDockerRunCommand(backendCall: BackendCall, image: String): String = {
    val callPath = containerCallPath(backendCall.workflowDescriptor, backendCall.call.unqualifiedName, backendCall.key.index)
    s"docker run --rm -v ${backendCall.callRootPath.toAbsolutePath}:$callPath -i $image"
  }


  private def runSubprocess(backendCall: BackendCall)(implicit ec: ExecutionContext): Future[ExecutionResult] = {
    val logger = workflowLoggerWithCall(backendCall)
    val stdoutWriter = backendCall.stdout.untailed
    val stderrTailed = backendCall.stderr.tailed(100)
    val dockerRun = backendCall.call.docker.map(d => buildDockerRunCommand(backendCall, d)).getOrElse("")
    val argv = Seq("/bin/bash", "-c", s"cat ${backendCall.script} | $dockerRun /bin/bash <&0")
    val process = argv.run(ProcessLogger(stdoutWriter writeWithNewline, stderrTailed writeWithNewline))

    // TODO: As currently implemented, this process.destroy() will kill the bash process but *not* its descendants. See ticket DSDEEPB-848.
    backendCall.callAbortRegistrationFunction.register(AbortFunction(() => process.destroy()))
    val backendCommandString = argv.map(s => "\""+s+"\"").mkString(" ")
    logger.info(s"command: $backendCommandString")
    val processReturnCode = process.exitValue() // blocks until process finishes
    Vector(stdoutWriter.writer, stderrTailed.writer) foreach { _.flushAndClose() }

    val stderrFileLength = Try(Files.size(backendCall.stderr)).getOrElse(0L)
    val returnCode = Try(
      if (processReturnCode == 0 || backendCall.call.docker.isEmpty) {
        val rc = backendCall.returnCode.contentAsString.stripLineEnd.toInt
        logger.info(s"Return code: $rc")
        rc
      } else {
        logger.error(s"Non-zero return code: $processReturnCode")
        logger.error(s"Standard error was:\n\n${stderrTailed.tailString}\n")
        throw new Exception(s"Unexpected process exit code: $processReturnCode")
      }
    )
    if (backendCall.call.failOnStderr && stderrFileLength > 0) {
      FailedExecution(new Throwable(s"Call ${backendCall.call.fullyQualifiedName}, " +
        s"Workflow ${backendCall.workflowDescriptor.id}: stderr has length $stderrFileLength")).future
    } else {

      def processSuccess(rc: Int) = {
        postProcess(backendCall) match {
          case Success(outputs) => backendCall.hash map { h => SuccessfulExecution(outputs, Seq.empty, rc, h) }
          case Failure(e) =>
            val message = Option(e.getMessage) map { ": " + _ } getOrElse ""
            FailedExecution(new Throwable("Failed post processing of outputs" + message, e)).future
        }
      }

      val badReturnCodeMessage =
         s"""Call ${backendCall.call.fullyQualifiedName}, Workflow ${backendCall.workflowDescriptor.id}: return code was ${returnCode.getOrElse("(none)")}
            |
            |Full command was: $backendCommandString
            |
            |${stderrTailed.tailString}
          """.stripMargin

      val continueOnReturnCode = backendCall.call.continueOnReturnCode
      returnCode match {
        case Success(143) => AbortedExecution.future // Special case to check for SIGTERM exit code - implying abort
        case Success(otherReturnCode) if continueOnReturnCode.continueFor(otherReturnCode) => processSuccess(otherReturnCode)
        case Success(badReturnCode) => FailedExecution(new Exception(badReturnCodeMessage), Option(badReturnCode)).future
        case Failure(e) => FailedExecution(new Throwable(badReturnCodeMessage, e)).future
      }
    }
  }
}
