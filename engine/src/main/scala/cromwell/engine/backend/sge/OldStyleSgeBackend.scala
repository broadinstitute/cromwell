package cromwell.engine.backend.sge

import java.nio.file.{Files, Path}

import akka.actor.ActorSystem
import better.files._
import com.google.api.client.util.ExponentialBackOff.Builder
import cromwell.CallEngineFunctions
import cromwell.backend.JobKey
import cromwell.engine._
import cromwell.engine.backend._
import cromwell.engine.backend.local.{OldStyleLocalBackend, SharedFileSystemBackend}
import cromwell.engine.backend.sge.OldStyleSgeBackend.InfoKeys
import cromwell.engine.db.DataAccess._
import cromwell.engine.workflow.BackendCallKey
import cromwell.logging.WorkflowLogger
import cromwell.util.FileUtil._

import scala.annotation.tailrec
import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.sys.process._
import scala.util.{Failure, Success, Try}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object OldStyleSgeBackend {
  object InfoKeys {
    val JobNumber = "SGE_JOB_NUMBER"
  }

  implicit class SgeEnhancedJobDescriptor(val jobDescriptor: OldStyleBackendCallJobDescriptor) extends AnyVal {
    def workflowRootPath = jobDescriptor.workflowDescriptor.workflowRootPath
    def stdout = jobDescriptor.callRootPath.resolve("stdout")
    def stderr = jobDescriptor.callRootPath.resolve("stderr")
    def script = jobDescriptor.callRootPath.resolve("script.sh")
    def returnCode = jobDescriptor.callRootPath.resolve("rc")
  }
}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class OldStyleSgeBackend(backendConfigEntry: BackendConfigurationEntry, actorSystem: ActorSystem) extends OldStyleBackend with SharedFileSystemBackend {
  def returnCode(jobDescriptor: OldStyleBackendCallJobDescriptor) = jobDescriptor.returnCode

  import OldStyleLocalBackend.WriteWithNewline

  override def backendType = BackendType.SGE

  override def adjustInputPaths(jobDescriptor: OldStyleBackendCallJobDescriptor) = adjustSharedInputPaths(jobDescriptor)

  /**
    * Exponential Backoff Builder to be used when polling for job status.
    */
  final private lazy val pollBackoffBuilder = new Builder()
    .setInitialIntervalMillis(10.seconds.toMillis.toInt)
    .setMaxElapsedTimeMillis(Int.MaxValue)
    .setMaxIntervalMillis(10.minutes.toMillis.toInt)
    .setMultiplier(1.1D)

  override def pollBackoff = pollBackoffBuilder.build()

  def stdoutStderr(jobDescriptor: OldStyleBackendCallJobDescriptor): CallLogs = sharedFileSystemStdoutStderr(jobDescriptor)

  def execute(jobDescriptor: OldStyleBackendCallJobDescriptor)(implicit ec: ExecutionContext): Future[OldStyleExecutionHandle] = Future( {
    val logger = jobLogger(jobDescriptor)
    jobDescriptor.instantiateCommand match {
      case Success(instantiatedCommand) =>
        logger.info(s"`$instantiatedCommand`")
        writeScript(jobDescriptor, instantiatedCommand)
        launchQsub(jobDescriptor) match {
          case (qsubReturnCode, _) if qsubReturnCode != 0 => OldStyleNonRetryableFailedExecution(
            new Throwable(s"Error: qsub exited with return code: $qsubReturnCode")).future
          case (_, None) => OldStyleNonRetryableFailedExecution(new Throwable(s"Could not parse Job ID from qsub output")).future
          case (_, Some(sgeJobId)) =>
            val updateDatabaseWithRunningInfo = updateSgeJobTable(jobDescriptor, "Running", None, Option(sgeJobId))
            pollForSgeJobCompletionThenPostProcess(jobDescriptor, sgeJobId) map { case (executionResult, jobRc) =>
              // Only send the completion update once the 'Running' update has completed (regardless of the success of that update)
              updateDatabaseWithRunningInfo onComplete {
                case _ =>
                  val completionStatus = statusString(executionResult)
                  val updateDatabaseWithCompletionInfo = updateSgeJobTable(jobDescriptor, completionStatus, Option(jobRc), Option(sgeJobId))
                  updateDatabaseWithCompletionInfo onFailure recordDatabaseFailure(logger, completionStatus, jobRc)
              }
              executionResult
            }
        }
      case Failure(ex) => OldStyleNonRetryableFailedExecution(ex).future
    }
  }).flatten map CompletedExecutionHandle

  private def statusString(result: OldStyleExecutionResult): String = (result match {
    case OldStyleAbortedExecution => ExecutionStatus.Aborted
    case OldStyleNonRetryableFailedExecution(_, _, _) => ExecutionStatus.Failed
    case OldStyleRetryableFailedExecution(_, _, _) => ExecutionStatus.Failed
    case OldStyleSuccessfulBackendCallExecution(_, _, _, _, _) => ExecutionStatus.Done
    case OldStyleSuccessfulFinalCallExecution => ExecutionStatus.Done
  }).toString

  private def recordDatabaseFailure(logger: WorkflowLogger, status: String, rc: Int): PartialFunction[Throwable, Unit] = {
    case e: Throwable =>
      logger.error(s"Failed to update database status: $status, rc: $rc because $e")
  }

  private def updateSgeJobTable(jobDescriptor: OldStyleBackendCallJobDescriptor, status: String, rc: Option[Int], sgeJobId: Option[Int])
                               (implicit ec: ExecutionContext): Future[Unit] = {
    globalDataAccess.updateExecutionInfo(jobDescriptor.workflowDescriptor.id, BackendCallKey(jobDescriptor.call, jobDescriptor.key.index, jobDescriptor.key.attempt), InfoKeys.JobNumber, Option(sgeJobId.toString))
  }

  /**
    * Returns the RC of this job when it finishes.  Sleeps and polls
    * until the 'rc' file is generated
    */
  private def waitUntilComplete(jobDescriptor: OldStyleBackendCallJobDescriptor): Int = {
    @tailrec
    def recursiveWait(): Int = Files.exists(jobDescriptor.returnCode) match {
      case true => jobDescriptor.returnCode.contentAsString.stripLineEnd.toInt
      case false =>
        val logger = jobLogger(jobDescriptor)
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
  private def killSgeJob(jobDescriptor: OldStyleBackendCallJobDescriptor, sgeJobId: Int) = () => {
    val qdelStdoutWriter = jobDescriptor.callRootPath.resolve("qdel.stdout").newBufferedWriter
    val qdelStderrWriter = jobDescriptor.callRootPath.resolve("qdel.stderr").newBufferedWriter
    val argv = Seq("qdel", sgeJobId.toString)
    val process = argv.run(ProcessLogger(qdelStdoutWriter writeWithNewline, qdelStderrWriter writeWithNewline))
    val returnCode: Int = process.exitValue()
    Vector(qdelStdoutWriter, qdelStderrWriter) foreach { _.flushAndClose() }
    jobDescriptor.returnCode.clear().appendLine("143")
    val logger = jobLogger(jobDescriptor)
    logger.debug(s"qdel $sgeJobId (returnCode=$returnCode)")
  }

  /**
    * Writes the script file containing the user's command from the WDL as well
    * as some extra shell code for monitoring jobs
    */
  private def writeScript(jobDescriptor: OldStyleBackendCallJobDescriptor, instantiatedCommand: String) = {
    jobDescriptor.script.write(
      s"""#!/bin/sh
          |$instantiatedCommand
          |echo $$? > rc
          |""".stripMargin)
  }

  /**
    * Launches the qsub command, returns a tuple: (rc, Option(sge_job_id))
    */
  private def launchQsub(jobDescriptor: OldStyleBackendCallJobDescriptor): (Int, Option[Int]) = {
    val logger = jobLogger(jobDescriptor)
    val sgeJobName = s"cromwell_${jobDescriptor.workflowDescriptor.shortId}_${jobDescriptor.call.unqualifiedName}"
    val argv = Seq("qsub", "-terse", "-N", sgeJobName, "-V", "-b", "n", "-wd", jobDescriptor.callRootPath.toAbsolutePath, "-o", jobDescriptor.stdout.getFileName, "-e", jobDescriptor.stderr.getFileName, jobDescriptor.script.toAbsolutePath).map(_.toString)
    val backendCommandString = argv.map(s => "\""+s+"\"").mkString(" ")
    logger.info(s"backend command: $backendCommandString")

    val qsubStdoutFile = jobDescriptor.callRootPath.resolve("qsub.stdout")
    val qsubStderrFile = jobDescriptor.callRootPath.resolve("qsub.stderr")
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
    * and returns the outputs for the jobDescriptor
    */
  private def pollForSgeJobCompletionThenPostProcess(jobDescriptor: OldStyleBackendCallJobDescriptor, sgeJobId: Int)(implicit ec: ExecutionContext): Future[(OldStyleExecutionResult, Int)] = {
    val logger = jobLogger(jobDescriptor)
    val abortFunction = killSgeJob(jobDescriptor, sgeJobId)
    val waitUntilCompleteFunction = waitUntilComplete(jobDescriptor)
    jobDescriptor.abortRegistrationFunction.foreach(_.register(AbortFunction(abortFunction)))
    val jobReturnCode = waitUntilCompleteFunction
    val continueOnReturnCode = jobDescriptor.callRuntimeAttributes.continueOnReturnCode
    logger.info(s"SGE job completed (returnCode=$jobReturnCode)")
    val executionResult = (jobReturnCode, jobDescriptor.stderr.toFile.length) match {
      case (r, _) if r == 143 => OldStyleAbortedExecution.future // Special case to check for SIGTERM exit code - implying abort
      case (r, _) if !continueOnReturnCode.continueFor(r) =>
        val message = s"SGE job failed because of return code: $r"
        logger.error(message)
        OldStyleNonRetryableFailedExecution(new Exception(message), Option(r)).future
      case (_, stderrLength) if stderrLength > 0 && jobDescriptor.callRuntimeAttributes.failOnStderr =>
        val message = s"SGE job failed because there were $stderrLength bytes on standard error"
        logger.error(message)
        OldStyleNonRetryableFailedExecution(new Exception(message), Option(0)).future
      case (r, _) =>
        postProcess(jobDescriptor) match {
          case Success(callOutputs) => jobDescriptor.hash map { h => OldStyleSuccessfulBackendCallExecution(callOutputs, Seq.empty, r, h) }
          case Failure(e) => OldStyleNonRetryableFailedExecution(e).future
        }
    }
    executionResult map { (_,  jobReturnCode) }
  }

  override def callRootPathWithBaseRoot(descriptor: OldStyleBackendCallJobDescriptor, baseRoot: String): Path = {
    val path = super.callRootPathWithBaseRoot(descriptor, baseRoot)
    if (!path.toFile.exists()) path.toFile.mkdirs()
    path
  }

  override def executionInfoKeys: List[String] = List(InfoKeys.JobNumber)

  override def callEngineFunctions(descriptor: OldStyleBackendCallJobDescriptor): CallEngineFunctions = {
    new SgeCallEngineFunctions(descriptor.workflowDescriptor.fileSystems, buildCallContext(descriptor))
  }

  def instantiateCommand(jobDescriptor: OldStyleBackendCallJobDescriptor): Try[String] = {
    val backendInputs = adjustInputPaths(jobDescriptor)
    jobDescriptor.call.instantiateCommandLine(backendInputs, jobDescriptor.callEngineFunctions)
  }

  override def poll(jobDescriptor: OldStyleBackendCallJobDescriptor, previous: OldStyleExecutionHandle)(implicit ec: ExecutionContext) = Future.successful(previous)

  override def resume(descriptor: OldStyleBackendCallJobDescriptor, executionInfos: Map[String, Option[String]])(implicit ec: ExecutionContext): Future[OldStyleExecutionHandle] = {
    Future.failed(new Throwable("resume invoked on non-resumable SGE backend"))
  }
  override def isResumable(key: JobKey, executionInfo: Map[String, Option[String]]) = false
  override def isRestartable(key: JobKey, executionInfo: Map[String, Option[String]]) = true
}
