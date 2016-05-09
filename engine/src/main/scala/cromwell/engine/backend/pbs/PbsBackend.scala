package cromwell.engine.backend.pbs

import java.nio.file.{Path, Files}

import akka.actor.ActorSystem
import better.files._
import com.google.api.client.util.ExponentialBackOff.Builder
import cromwell.engine.backend._
import cromwell.engine.backend.local.{LocalBackend, SharedFileSystem}
import cromwell.engine.backend.pbs.PbsBackend.InfoKeys
import cromwell.engine.db.DataAccess._
import cromwell.engine.workflow.BackendCallKey
import cromwell.engine.{AbortRegistrationFunction, _}
import cromwell.logging.WorkflowLogger
import cromwell.util.FileUtil._
import cromwell.util.TryUtil

import scala.annotation.tailrec
import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.sys.process._
import scala.util.{Failure, Success, Try}

object PbsBackend {
  object InfoKeys {
    val JobNumber = "PBS_JOB_NUMBER"
  }
}

case class PbsBackend(actorSystem: ActorSystem) extends Backend with SharedFileSystem {
  type BackendCall = PbsBackendCall

  import LocalBackend.WriteWithNewline

  override def backendType = BackendType.PBS

  override def adjustInputPaths(jobDescriptor: BackendCallJobDescriptor) = adjustSharedInputPaths(jobDescriptor)

  /**
    * Exponential Backoff Builder to be used when polling for call status.
    */
  final private lazy val pollBackoffBuilder = new Builder()
    .setInitialIntervalMillis(10.seconds.toMillis.toInt)
    .setMaxElapsedTimeMillis(Int.MaxValue)
    .setMaxIntervalMillis(10.minutes.toMillis.toInt)
    .setMultiplier(1.1D)

  override def pollBackoff = pollBackoffBuilder.build()

  override def bindCall(jobDescriptor: BackendCallJobDescriptor,
                        abortRegistrationFunction: Option[AbortRegistrationFunction]): BackendCall = {
    PbsBackendCall(this, jobDescriptor, abortRegistrationFunction)
  }

  def stdoutStderr(backendCall: BackendCall): CallLogs = sharedFileSystemStdoutStderr(backendCall)

  def execute(backendCall: BackendCall)(implicit ec: ExecutionContext): Future[ExecutionHandle] = Future( {
    val logger = workflowLoggerWithCall(backendCall)
    backendCall.instantiateCommand match {
      case Success(instantiatedCommand) =>
        logger.info(s"`$instantiatedCommand`")
        writeScript(backendCall, instantiatedCommand)
        launchQsub(backendCall) match {
          case (qsubReturnCode, _) if qsubReturnCode != 0 => NonRetryableExecution(
            new Throwable(s"Error: qsub exited with return code: $qsubReturnCode")).future
          case (_, None) => NonRetryableExecution(new Throwable(s"Could not parse Job ID from qsub output")).future
          case (_, Some(pbsJobId)) =>
            val updateDatabaseWithRunningInfo = updatePbsJobTable(backendCall, "Running", None, Option(pbsJobId))
            pollForPbsJobCompletionThenPostProcess(backendCall, pbsJobId) map { case (executionResult, jobRc) =>
              // Only send the completion update once the 'Running' update has completed (regardless of the success of that update)
              updateDatabaseWithRunningInfo onComplete {
                case _ =>
                  val completionStatus = statusString(executionResult)
                  val updateDatabaseWithCompletionInfo = updatePbsJobTable(backendCall, completionStatus, Option(jobRc), Option(pbsJobId))
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

  private def updatePbsJobTable(call: BackendCall, status: String, rc: Option[Int], pbsJobId: Option[Int])
                               (implicit ec: ExecutionContext): Future[Unit] = {
    globalDataAccess.updateExecutionInfo(call.workflowDescriptor.id, BackendCallKey(call.call, call.key.index, call.key.attempt), InfoKeys.JobNumber, Option(pbsJobId.toString))
  }

  /** TODO restart isn't currently implemented for PBS, there is probably work that needs to be done here much like
    * JES restart, which perhaps could be factored out into a common "remote executor" trait.
    */
  override def prepareForRestart(restartableWorkflow: WorkflowDescriptor)(implicit ec: ExecutionContext) = Future.successful(())

  /**
   * Returns the RC of this job when it finishes.  Sleeps and polls until the 'rc' file is
   * generated and stdout, stderr files exist in their specified locations: a feature/quirk
   * of PBS is that it creates output stream files in /var/spool/PBS/spool/... and only after
   * execution completes are they are copied to the locations specified by -o/-e.
   */
  private def waitUntilComplete(backendCall: BackendCall): Int = {
    @tailrec
    def recursiveWait(): Int = Seq(backendCall.returnCode, backendCall.stdout, backendCall.stderr).forall(Files.exists(_)) match {
      case true => backendCall.returnCode.contentAsString.stripLineEnd.toInt
      case false =>
        val logger = workflowLoggerWithCall(backendCall)
        logger.info("output files do not all exist yet")
        Thread.sleep(5000)
        recursiveWait()
    }
    recursiveWait()
  }

  /**
   * This returns a zero-parameter function which, when called, will kill the
   * PBS job with id `pbsJobId`.  It also writes to the 'rc' file the value '143'
   */
  private def killPbsJob(backendCall: BackendCall, pbsJobId: Int) = () => {
    val qdelStdoutWriter = backendCall.callRootPath.resolve("qdel.stdout").newBufferedWriter
    val qdelStderrWriter = backendCall.callRootPath.resolve("qdel.stderr").newBufferedWriter
    val argv = Seq("qdel", pbsJobId.toString)
    val process = argv.run(ProcessLogger(qdelStdoutWriter writeWithNewline, qdelStderrWriter writeWithNewline))
    val returnCode: Int = process.exitValue()
    Vector(qdelStdoutWriter, qdelStderrWriter) foreach { _.flushAndClose() }
    backendCall.returnCode.clear().appendLine("143")
    val logger = workflowLoggerWithCall(backendCall)
    logger.debug(s"qdel $pbsJobId (returnCode=$returnCode)")
  }

  /**
   * Writes the script file containing the user's command from the WDL as well
   * as some extra shell code for monitoring jobs
   */
  private def writeScript(backendCall: BackendCall, instantiatedCommand: String) = {
    backendCall.script.write(
      /*
       * Use ':' as stripMargin token, as we've had pipe '|' appear legitimately
       * at the beginning of a backslash-continued line of shell script in
       * $instantiatedCommand and get stripped out. We don't want that! Colon
       * seems a kind of unlikely thing to appear there, additionally, ':' in
       * bash is a no-op so it's probably safe to remove on the rare occasions
       * it does occur at start of a line. We could just forget about it and
       * have an ugly left margin-abutted embedded string but this seems a nice
       * compromise even if it has required a long comment to justify... haha.
       */
      s"""#!/bin/sh
         :cd $$PBS_O_WORKDIR
         :$instantiatedCommand
         :echo $$? > rc
         :""".stripMargin(':'))
  }

  /**
   * Launches the qsub command, returns a tuple: (rc, Option(pbs_job_id))
   */
  private def launchQsub(backendCall: BackendCall): (Int, Option[Int]) = {
    val logger = workflowLoggerWithCall(backendCall)
    val pbsJobName = s"${backendCall.workflowDescriptor.shortId}_${backendCall.call.unqualifiedName}" take 15
    val queueSpec = backendCall.runtimeAttributes.queue match {
      case Some(q) => List("-q", q)
      case None => List()
    }
    val argv = "qsub" :: queueSpec ++ Seq(
      "-N", pbsJobName,
      "-e", backendCall.stderr.toAbsolutePath,
      "-o", backendCall.stdout.toAbsolutePath,
      "-l", s"ncpus=${backendCall.runtimeAttributes.cpu}",
      "-l", s"mem=${backendCall.runtimeAttributes.memoryGB.toInt}gb",
      "-l", s"walltime=${backendCall.runtimeAttributes.walltime}",
      backendCall.script.toAbsolutePath
    ).map(_.toString)
    val backendCommandString = argv.map(s => "\""+s+"\"").mkString(" ")
    logger.info(s"backend command: $backendCommandString")

    val qsubStdoutFile = backendCall.callRootPath.resolve("qsub.stdout")
    val qsubStderrFile = backendCall.callRootPath.resolve("qsub.stderr")
    val qsubStdoutWriter = qsubStdoutFile.newBufferedWriter
    val qsubStderrWriter = qsubStderrFile.newBufferedWriter
    val executionDir = backendCall.callRootPath.toAbsolutePath.toFile
    val process = Process(argv, executionDir).run(ProcessLogger(qsubStdoutWriter writeWithNewline, qsubStderrWriter writeWithNewline))
    val returnCode: Int = process.exitValue()
    Vector(qsubStdoutWriter, qsubStderrWriter) foreach { _.flushAndClose() }

    // stdout only has the job ID, if it was successfully launched
    val jobId = Try(qsubStdoutFile.contentAsString.stripLineEnd.split("\\.")(0).toInt) match {
      case Success(id) => Some(id)
      case Failure(ex) =>
        logger.error(s"Could not find PBS job ID from qsub stdout file.\n\n" +
          s"Check the qsub stderr file for possible errors: ${qsubStderrFile.toAbsolutePath}")
        None
    }
    logger.info(s"qsub returnCode=$returnCode, job ID=${jobId.getOrElse("NONE")}")
    (returnCode, jobId)
  }

  /**
   * This waits for a given PBS job to finish.  When finished, it post-processes the job
   * and returns the outputs for the call
   */
  private def pollForPbsJobCompletionThenPostProcess(backendCall: BackendCall, pbsJobId: Int)(implicit ec: ExecutionContext): Future[(ExecutionResult, Int)] = {
    val logger = workflowLoggerWithCall(backendCall)
    val abortFunction = killPbsJob(backendCall, pbsJobId)
    val waitUntilCompleteFunction = waitUntilComplete(backendCall)
    backendCall.callAbortRegistrationFunction.foreach(_.register(AbortFunction(abortFunction)))
    val jobReturnCode = waitUntilCompleteFunction
    val continueOnReturnCode = backendCall.runtimeAttributes.continueOnReturnCode
    logger.info(s"PBS job completed (returnCode=$jobReturnCode)")
    val executionResult = (jobReturnCode, backendCall.stderr.toFile.length) match {
      case (r, _) if r == 143 => AbortedExecution.future // Special case to check for SIGTERM exit code - implying abort
      case (r, _) if !continueOnReturnCode.continueFor(r) =>
        val message = s"PBS job failed because of return code: $r"
        logger.error(message)
        NonRetryableExecution(new Exception(message), Option(r)).future
      case (_, stderrLength) if stderrLength > 0 && backendCall.runtimeAttributes.failOnStderr =>
        val message = s"PBS job failed because there were $stderrLength bytes on standard error"
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

  override def callRootPathWithBaseRoot(jobDescriptor: BackendCallJobDescriptor, baseRoot: String): Path = {
    val path = super.callRootPathWithBaseRoot(jobDescriptor, baseRoot)
    if (!path.toFile.exists()) path.toFile.mkdirs()
    path
  }

  override def executionInfoKeys: List[String] = List(InfoKeys.JobNumber)

  override def callEngineFunctions(descriptor: BackendCallJobDescriptor): CallEngineFunctions = {
    new PbsCallEngineFunctions(descriptor.workflowDescriptor.ioManager, buildCallContext(descriptor))
  }

  def instantiateCommand(jobDescriptor: BackendCallJobDescriptor): Try[String] = {
    val backendInputs = adjustInputPaths(jobDescriptor)
    jobDescriptor.key.scope.instantiateCommandLine(backendInputs, jobDescriptor.callEngineFunctions)
  }

  override def poll(jobDescriptor: BackendCallJobDescriptor, previous: ExecutionHandle)(implicit ec: ExecutionContext) = Future.successful(previous)
}
