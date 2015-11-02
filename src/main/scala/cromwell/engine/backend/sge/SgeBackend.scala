package cromwell.engine.backend.sge

import java.nio.file.Files

import com.typesafe.scalalogging.LazyLogging
import cromwell.binding.{Call, CallInputs}
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.backend._
import cromwell.engine.backend.local.{LocalBackend, SharedFileSystem}
import cromwell.engine.db.DataAccess._
import cromwell.engine.db.SgeCallBackendInfo
import cromwell.engine.workflow.CallKey
import cromwell.engine.{AbortRegistrationFunction, _}
import cromwell.parser.BackendType
import cromwell.util.FileUtil._

import scala.annotation.tailrec
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.sys.process._
import scala.util.{Failure, Success, Try}

class SgeBackend extends Backend with SharedFileSystem with LazyLogging {
  type BackendCall = SgeBackendCall

  import LocalBackend.WriteWithNewline

  override def backendType = BackendType.SGE

  override def bindCall(workflowDescriptor: WorkflowDescriptor,
                        key: CallKey,
                        locallyQualifiedInputs: CallInputs,
                        abortRegistrationFunction: AbortRegistrationFunction): BackendCall = {
    SgeBackendCall(this, workflowDescriptor, key, locallyQualifiedInputs, abortRegistrationFunction)
  }

  def execute(backendCall: BackendCall)(implicit ec: ExecutionContext): Future[ExecutionHandle] = Future {
    val tag = makeTag(backendCall)
    backendCall.instantiateCommand match {
      case Success(instantiatedCommand) =>
        logger.info(s"$tag `$instantiatedCommand`")
        writeScript(backendCall, instantiatedCommand)
        launchQsub(backendCall) match {
          case (qsubReturnCode, _) if qsubReturnCode != 0 => FailedExecution(
            new Throwable(s"Error: qsub exited with return code: $qsubReturnCode"))
          case (_, None) => FailedExecution(new Throwable(s"Could not parse Job ID from qsub output"))
          case (_, Some(sgeJobId)) =>
            val updateDatabaseWithRunningInfo = updateSgeJobTable(backendCall, "Running", None, Option(sgeJobId))
            val (executionResult, jobRc) = pollForSgeJobCompletionThenPostProcess(backendCall, sgeJobId)
            // Only send the completion update once the 'Running' update has completed (regardless of the success of that update)
            updateDatabaseWithRunningInfo onComplete {
              case _ =>
                val completionStatus = statusString(executionResult)
                val updateDatabaseWithCompletionInfo = updateSgeJobTable(backendCall, completionStatus, Option(jobRc), Option(sgeJobId))
                updateDatabaseWithCompletionInfo onFailure recordDatabaseFailure(backendCall.call, completionStatus, jobRc)
            }

            executionResult
        }
      case Failure(ex) => FailedExecution(ex)
    }
  } map CompletedExecutionHandle

  private def statusString(result: ExecutionResult): String = (result match {
      case AbortedExecution => ExecutionStatus.Aborted
      case FailedExecution(_, _) => ExecutionStatus.Failed
      case SuccessfulExecution(_, _) => ExecutionStatus.Done
    }).toString

  private def recordDatabaseFailure(call: Call, status: String, rc: Int): PartialFunction[Throwable, Unit] = {
    case e: Throwable => logger.error(s"Failed to update database with ${call.name} status: $status, rc: $rc because $e")
  }

  private def updateSgeJobTable(call: BackendCall, status: String, rc: Option[Int], sgeJobId: Option[Int]): Future[Unit] = {
    val backendInfo = SgeCallBackendInfo(sgeJobId)
    globalDataAccess.updateExecutionBackendInfo(call.workflowDescriptor.id, CallKey(call.call, call.key.index), backendInfo)
  }

  /** TODO restart isn't currently implemented for SGE, there is probably work that needs to be done here much like
    * JES restart, which perhaps could be factored out into a common "remote executor" trait.
    */
  override def prepareForRestart(restartableWorkflow: WorkflowDescriptor)(implicit ec: ExecutionContext) = Future.successful(())

  /**
   * Returns the RC of this job when it finishes.  Sleeps and polls
   * until the 'rc' file is generated
   */
  private def waitUntilComplete(backendCall: BackendCall): Int = {
    val tag = makeTag(backendCall)
    @tailrec
    def recursiveWait(): Int = Files.exists(backendCall.returnCode) match {
      case true => backendCall.returnCode.toFile.slurp.stripLineEnd.toInt
      case false =>
        logger.info(s"$tag 'rc' file does not exist yet")
        Thread.sleep(5000)
        recursiveWait()
    }
    recursiveWait()
  }

  /**
   * This returns a zero-parameter function which, when called, will kill the
   * SGE job with id `sgeJobId`.  It also writes to the 'rc' file the value '143'
   */
  private def killSgeJob(backendCall: BackendCall, sgeJobId: Int) = () => {
    val (_, qdelStdoutWriter) = backendCall.callRootPath.resolve("qdel.stdout").fileAndWriter
    val (_, qdelStderrWriter) = backendCall.callRootPath.resolve("qdel.stderr").fileAndWriter
    val argv = Seq("qdel", sgeJobId.toString)
    val process = argv.run(ProcessLogger(qdelStdoutWriter writeWithNewline, qdelStderrWriter writeWithNewline))
    val returnCode: Int = process.exitValue()
    Vector(qdelStdoutWriter, qdelStderrWriter) foreach { _.flushAndClose() }
    val (_, rcWriter) = backendCall.returnCode.fileAndWriter
    rcWriter.writeWithNewline("143")
    rcWriter.flushAndClose()
    logger.debug(s"qdel $sgeJobId (returnCode=$returnCode)")
  }

  /**
   * Writes the script file containing the user's command from the WDL as well
   * as some extra shell code for monitoring jobs
   */
  private def writeScript(backendCall: BackendCall, instantiatedCommand: String) = {
    val (_, scriptWriter) = backendCall.script.fileAndWriter
    scriptWriter.writeWithNewline("#!/bin/sh")
    scriptWriter.writeWithNewline(instantiatedCommand)
    scriptWriter.writeWithNewline("echo $? > rc")
    scriptWriter.flushAndClose()
  }

  /**
   * Launches the qsub command, returns a tuple: (rc, Option(sge_job_id))
   */
  private def launchQsub(backendCall: BackendCall): (Int, Option[Int]) = {
    val tag = makeTag(backendCall)
    val sgeJobName = s"cromwell_${backendCall.workflowDescriptor.shortId}_${backendCall.call.name}"
    val argv = Seq("qsub", "-terse", "-N", sgeJobName, "-V", "-b", "n", "-wd", backendCall.callRootPath.toAbsolutePath, "-o", backendCall.stdout.getFileName, "-e", backendCall.stderr.getFileName, backendCall.script.toAbsolutePath).map(_.toString)
    val backendCommandString = argv.map(s => "\""+s+"\"").mkString(" ")
    logger.info(s"$tag backend command: $backendCommandString")

    val (qsubStdoutFile, qsubStdoutWriter) = backendCall.callRootPath.resolve("qsub.stdout").fileAndWriter
    val (qsubStderrFile, qsubStderrWriter) = backendCall.callRootPath.resolve("qsub.stderr").fileAndWriter
    val process = argv.run(ProcessLogger(qsubStdoutWriter writeWithNewline, qsubStderrWriter writeWithNewline))
    val returnCode: Int = process.exitValue()
    Vector(qsubStdoutWriter, qsubStderrWriter) foreach { _.flushAndClose() }

    // The -terse option to qsub makes it so stdout only has the job ID, if it was successfully launched
    val jobId = Try(qsubStdoutFile.slurp.stripLineEnd.toInt) match {
      case Success(id) => Some(id)
      case Failure(ex) =>
        logger.error(s"$tag Could not find SGE job ID from qsub stdout file.\n\nCheck the qsub stderr file for possible errors: ${qsubStderrFile.toAbsolutePath.toString}")
        None
    }
    logger.info(s"$tag qsub returnCode=$returnCode, job ID=${jobId.getOrElse("NONE")}")
    (returnCode, jobId)
  }

  /**
   * This waits for a given SGE job to finish.  When finished, it post-processes the job
   * and returns the outputs for the call
   */
  private def pollForSgeJobCompletionThenPostProcess(backendCall: BackendCall, sgeJobId: Int): (ExecutionResult, Int) = {
    val tag = makeTag(backendCall)
    val abortFunction = killSgeJob(backendCall, sgeJobId)
    val waitUntilCompleteFunction = waitUntilComplete(backendCall)
    backendCall.callAbortRegistrationFunction.register(AbortFunction(abortFunction))
    val jobReturnCode = waitUntilCompleteFunction
    val continueOnReturnCode = backendCall.call.continueOnReturnCode
    logger.info(s"$tag SGE job completed (returnCode=$jobReturnCode)")
    val executionResult = (jobReturnCode, backendCall.stderr.toFile.length) match {
      case (r, _) if r == 143 => AbortedExecution // Special case to check for SIGTERM exit code - implying abort
      case (r, _) if !continueOnReturnCode.continueFor(r) =>
        FailedExecution(new Exception(s"$tag SGE job failed because of return code: $r"), Option(r))
      case (_, stderrLength) if stderrLength > 0 && backendCall.call.failOnStderr =>
        FailedExecution(new Throwable(s"$tag SGE job failed because there were $stderrLength bytes on standard error"), Option(0))
      case (r, _) =>
        postProcess(backendCall) match {
          case Success(callOutputs) => SuccessfulExecution(callOutputs, r)
          case Failure(e) => FailedExecution(e)
        }
    }
    (executionResult, jobReturnCode)
  }
}
