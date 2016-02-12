package cromwell.engine.backend.sge

import java.nio.file.Files

import akka.actor.ActorSystem
import better.files._
import com.google.api.client.util.ExponentialBackOff.Builder
import cromwell.engine.backend._
import cromwell.engine.backend.local.{LocalBackend, SharedFileSystem}
import cromwell.engine.backend.sge.SgeBackend.InfoKeys
import cromwell.engine.db.DataAccess._
import cromwell.engine.workflow.BackendCallKey
import cromwell.engine.{AbortRegistrationFunction, _}
import cromwell.logging.WorkflowLogger
import cromwell.util.FileUtil._
import wdl4s.CallInputs

import scala.annotation.tailrec
import scala.concurrent.{ExecutionContext, Future}
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.sys.process._
import scala.util.{Failure, Success, Try}

object SgeBackend {
  object InfoKeys {
    val JobNumber = "SGE_JOB_NUMBER"
  }
}

case class SgeBackend(actorSystem: ActorSystem) extends Backend with SharedFileSystem {
  type BackendCall = SgeBackendCall

  import LocalBackend.WriteWithNewline

  override def backendType = BackendType.SGE

  override def adjustInputPaths(backendCall: BackendCall) = adjustSharedInputPaths(backendCall)

  /**
    * Exponential Backoff Builder to be used when polling for call status.
    */
  final private lazy val pollBackoffBuilder = new Builder()
    .setInitialIntervalMillis(10.seconds.toMillis.toInt)
    .setMaxElapsedTimeMillis(Int.MaxValue)
    .setMaxIntervalMillis(10.minutes.toMillis.toInt)
    .setMultiplier(1.1D)

  override def pollBackoff = pollBackoffBuilder.build()

  override def bindCall(workflowDescriptor: WorkflowDescriptor,
                        key: BackendCallKey,
                        locallyQualifiedInputs: CallInputs,
                        abortRegistrationFunction: Option[AbortRegistrationFunction]): BackendCall = {
    SgeBackendCall(this, workflowDescriptor, key, locallyQualifiedInputs, abortRegistrationFunction)
  }

  def stdoutStderr(backendCall: BackendCall): CallLogs = sharedFileSystemStdoutStderr(backendCall)

  def execute(backendCall: BackendCall)(implicit ec: ExecutionContext): Future[ExecutionHandle] = Future( {
    val logger = workflowLoggerWithCall(backendCall)
    instantiateCommand(backendCall) match {
      case Success(instantiatedCommand) =>
        logger.info(s"`$instantiatedCommand`")
        writeScript(backendCall, instantiatedCommand)
        launchQsub(backendCall) match {
          case (qsubReturnCode, _) if qsubReturnCode != 0 => NonRetryableExecution(
            new Throwable(s"Error: qsub exited with return code: $qsubReturnCode")).future
          case (_, None) => NonRetryableExecution(new Throwable(s"Could not parse Job ID from qsub output")).future
          case (_, Some(sgeJobId)) =>
            val updateDatabaseWithRunningInfo = updateSgeJobTable(backendCall, "Running", None, Option(sgeJobId))
            pollForSgeJobCompletionThenPostProcess(backendCall, sgeJobId) map { case (executionResult, jobRc) =>
              // Only send the completion update once the 'Running' update has completed (regardless of the success of that update)
              updateDatabaseWithRunningInfo onComplete {
                case _ =>
                  val completionStatus = statusString(executionResult)
                  val updateDatabaseWithCompletionInfo = updateSgeJobTable(backendCall, completionStatus, Option(jobRc), Option(sgeJobId))
                  updateDatabaseWithCompletionInfo onFailure recordDatabaseFailure(logger, completionStatus, jobRc)
              }
              executionResult
            }
        }
      case Failure(ex) => NonRetryableExecution(ex).future
    }
  }).flatten map CompletedExecutionHandle

  private def statusString(result: ExecutionResult): String = (result match {
      case AbortedExecution => ExecutionStatus.Aborted
      case NonRetryableExecution(_, _, _) => ExecutionStatus.Failed
      case RetryableExecution(_, _, _) => ExecutionStatus.Failed
      case SuccessfulBackendCallExecution(_, _, _, _, _) => ExecutionStatus.Done
      case SuccessfulFinalCallExecution => ExecutionStatus.Done
    }).toString

  private def recordDatabaseFailure(logger: WorkflowLogger, status: String, rc: Int): PartialFunction[Throwable, Unit] = {
    case e: Throwable =>
      logger.error(s"Failed to update database status: $status, rc: $rc because $e")
  }

  private def updateSgeJobTable(call: BackendCall, status: String, rc: Option[Int], sgeJobId: Option[Int])
                               (implicit ec: ExecutionContext): Future[Unit] = {
    globalDataAccess.updateExecutionInfo(call.workflowDescriptor.id, BackendCallKey(call.call, call.key.index, call.key.attempt), InfoKeys.JobNumber, Option(sgeJobId.toString))
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
    @tailrec
    def recursiveWait(): Int = Files.exists(backendCall.returnCode) match {
      case true => backendCall.returnCode.contentAsString.stripLineEnd.toInt
      case false =>
        val logger = workflowLoggerWithCall(backendCall)
        logger.info(s"'rc' file does not exist yet")
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
    val qdelStdoutWriter = backendCall.callRootPath.resolve("qdel.stdout").newBufferedWriter
    val qdelStderrWriter = backendCall.callRootPath.resolve("qdel.stderr").newBufferedWriter
    val argv = Seq("qdel", sgeJobId.toString)
    val process = argv.run(ProcessLogger(qdelStdoutWriter writeWithNewline, qdelStderrWriter writeWithNewline))
    val returnCode: Int = process.exitValue()
    Vector(qdelStdoutWriter, qdelStderrWriter) foreach { _.flushAndClose() }
    backendCall.returnCode.clear().appendLine("143")
    val logger = workflowLoggerWithCall(backendCall)
    logger.debug(s"qdel $sgeJobId (returnCode=$returnCode)")
  }

  /**
   * Writes the script file containing the user's command from the WDL as well
   * as some extra shell code for monitoring jobs
   */
  private def writeScript(backendCall: BackendCall, instantiatedCommand: String) = {
    backendCall.script.write(
      s"""#!/bin/sh
         |$instantiatedCommand
         |echo $$? > rc
         |""".stripMargin)
  }

  /**
   * Launches the qsub command, returns a tuple: (rc, Option(sge_job_id))
   */
  private def launchQsub(backendCall: BackendCall): (Int, Option[Int]) = {
    val logger = workflowLoggerWithCall(backendCall)
    val sgeJobName = s"cromwell_${backendCall.workflowDescriptor.shortId}_${backendCall.call.unqualifiedName}"
    val argv = Seq("qsub", "-terse", "-N", sgeJobName, "-V", "-b", "n", "-wd", backendCall.callRootPath.toAbsolutePath, "-o", backendCall.stdout.getFileName, "-e", backendCall.stderr.getFileName, backendCall.script.toAbsolutePath).map(_.toString)
    val backendCommandString = argv.map(s => "\""+s+"\"").mkString(" ")
    logger.info(s"backend command: $backendCommandString")

    val qsubStdoutFile = backendCall.callRootPath.resolve("qsub.stdout")
    val qsubStderrFile = backendCall.callRootPath.resolve("qsub.stderr")
    val qsubStdoutWriter = qsubStdoutFile.newBufferedWriter
    val qsubStderrWriter = qsubStderrFile.newBufferedWriter
    val process = argv.run(ProcessLogger(qsubStdoutWriter writeWithNewline, qsubStderrWriter writeWithNewline))
    val returnCode: Int = process.exitValue()
    Vector(qsubStdoutWriter, qsubStderrWriter) foreach { _.flushAndClose() }

    // The -terse option to qsub makes it so stdout only has the job ID, if it was successfully launched
    val jobId = Try(qsubStdoutFile.contentAsString.stripLineEnd.toInt) match {
      case Success(id) => Some(id)
      case Failure(ex) =>
        logger.error(s"Could not find SGE job ID from qsub stdout file.\n\n" +
          s"Check the qsub stderr file for possible errors: ${qsubStderrFile.toAbsolutePath}")
        None
    }
    logger.info(s"qsub returnCode=$returnCode, job ID=${jobId.getOrElse("NONE")}")
    (returnCode, jobId)
  }

  /**
   * This waits for a given SGE job to finish.  When finished, it post-processes the job
   * and returns the outputs for the call
   */
  private def pollForSgeJobCompletionThenPostProcess(backendCall: BackendCall, sgeJobId: Int)(implicit ec: ExecutionContext): Future[(ExecutionResult, Int)] = {
    val logger = workflowLoggerWithCall(backendCall)
    val abortFunction = killSgeJob(backendCall, sgeJobId)
    val waitUntilCompleteFunction = waitUntilComplete(backendCall)
    backendCall.callAbortRegistrationFunction.foreach(_.register(AbortFunction(abortFunction)))
    val jobReturnCode = waitUntilCompleteFunction
    val continueOnReturnCode = backendCall.runtimeAttributes.continueOnReturnCode
    logger.info(s"SGE job completed (returnCode=$jobReturnCode)")
    val executionResult = (jobReturnCode, backendCall.stderr.toFile.length) match {
      case (r, _) if r == 143 => AbortedExecution.future // Special case to check for SIGTERM exit code - implying abort
      case (r, _) if !continueOnReturnCode.continueFor(r) =>
        val message = s"SGE job failed because of return code: $r"
        logger.error(message)
        NonRetryableExecution(new Exception(message), Option(r)).future
      case (_, stderrLength) if stderrLength > 0 && backendCall.runtimeAttributes.failOnStderr =>
        val message = s"SGE job failed because there were $stderrLength bytes on standard error"
        logger.error(message)
        NonRetryableExecution(new Exception(message), Option(0)).future
      case (r, _) =>
        postProcess(backendCall) match {
          case Success(callOutputs) => backendCall.hash map { h => SuccessfulBackendCallExecution(callOutputs, Seq.empty, r, h) }
          case Failure(e) => NonRetryableExecution(e).future
        }
    }
    executionResult map { (_,  jobReturnCode) }
  }

  override def executionInfoKeys: List[String] = List(InfoKeys.JobNumber)
}
