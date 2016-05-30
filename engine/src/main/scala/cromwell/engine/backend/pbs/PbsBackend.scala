package cromwell.engine.backend.pbs

import java.nio.file.{FileSystem, Files, Path}

import akka.actor.ActorSystem
import better.files._
import com.google.api.client.util.ExponentialBackOff.Builder
import cromwell.core.WorkflowOptions
import cromwell.engine._
import cromwell.engine.backend._
import cromwell.engine.backend.io._
import cromwell.engine.backend.io.filesystem.gcs.{GcsFileSystemProvider, StorageFactory}
import cromwell.engine.backend.local.{LocalBackend, SharedFileSystem}
import cromwell.engine.backend.pbs.PbsBackend.InfoKeys
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

object PbsBackend {
  object InfoKeys {
    val JobNumber = "PBS_JOB_NUMBER"
  }

  implicit class PbsEnhancedJobDescriptor(val jobDescriptor: BackendCallJobDescriptor) extends AnyVal {
    def workflowRootPath = jobDescriptor.workflowDescriptor.workflowRootPath
    def stdout = jobDescriptor.callRootPath.resolve("stdout")
    def stderr = jobDescriptor.callRootPath.resolve("stderr")
    def script = jobDescriptor.callRootPath.resolve("script.sh")
    def returnCode = jobDescriptor.callRootPath.resolve("rc")
  }
}

case class PbsBackend(actorSystem: ActorSystem) extends Backend with SharedFileSystem {
  val installationIsPBSPro: Boolean = ("qstat --version" #| "grep PBSPro" !) == 0
  
  def returnCode(jobDescriptor: BackendCallJobDescriptor) = jobDescriptor.returnCode

  import LocalBackend.WriteWithNewline

  override def backendType = BackendType.PBS

  override def adjustInputPaths(jobDescriptor: BackendCallJobDescriptor) = adjustSharedInputPaths(jobDescriptor)

  /**
    * Exponential Backoff Builder to be used when polling for job status.
    */
  final private lazy val pollBackoffBuilder = new Builder()
    .setInitialIntervalMillis(10.seconds.toMillis.toInt)
    .setMaxElapsedTimeMillis(Int.MaxValue)
    .setMaxIntervalMillis(10.minutes.toMillis.toInt)
    .setMultiplier(1.1D)

  override def pollBackoff = pollBackoffBuilder.build()

  def stdoutStderr(jobDescriptor: BackendCallJobDescriptor): CallLogs = sharedFileSystemStdoutStderr(jobDescriptor)

  def execute(jobDescriptor: BackendCallJobDescriptor)(implicit ec: ExecutionContext): Future[ExecutionHandle] = Future( {
    val logger = jobLogger(jobDescriptor)
    jobDescriptor.instantiateCommand match {
      case Success(instantiatedCommand) =>
        logger.info(s"`$instantiatedCommand`")
        writeScript(jobDescriptor, instantiatedCommand)
        launchQsub(jobDescriptor) match {
          case (qsubReturnCode, _) if qsubReturnCode != 0 => NonRetryableExecution(
            new Throwable(s"Error: qsub exited with return code: $qsubReturnCode")).future
          case (_, None) => NonRetryableExecution(new Throwable(s"Could not parse Job ID from qsub output")).future
          case (_, Some(pbsJobId)) =>
            val updateDatabaseWithRunningInfo = updatePbsJobTable(jobDescriptor, "Running", None, Option(pbsJobId))
            pollForPbsJobCompletionThenPostProcess(jobDescriptor, pbsJobId) map { case (executionResult, jobRc) =>
              // Only send the completion update once the 'Running' update has completed (regardless of the success of that update)
              updateDatabaseWithRunningInfo onComplete {
                case _ =>
                  val completionStatus = statusString(executionResult)
                  val updateDatabaseWithCompletionInfo = updatePbsJobTable(jobDescriptor, completionStatus, Option(jobRc), Option(pbsJobId))
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

  private def updatePbsJobTable(jobDescriptor: BackendCallJobDescriptor, status: String, rc: Option[Int], pbsJobId: Option[Int])
                               (implicit ec: ExecutionContext): Future[Unit] = {
    globalDataAccess.updateExecutionInfo(jobDescriptor.workflowDescriptor.id, BackendCallKey(jobDescriptor.call, jobDescriptor.key.index, jobDescriptor.key.attempt), InfoKeys.JobNumber, Option(pbsJobId.toString))
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
  private def waitUntilComplete(jobDescriptor: BackendCallJobDescriptor, pbsJobId: Int): Int = {
    val logger = jobLogger(jobDescriptor)
    @tailrec
    def recursiveWait(): Int = (Seq(jobDescriptor.returnCode, jobDescriptor.stdout, jobDescriptor.stderr).forall(Files.exists(_)), jobDone(pbsJobId)) match {
      case (true, _) => jobDescriptor.returnCode.contentAsString.stripLineEnd.toInt
      case (false, false) =>
        logger.info(s"Output files do not exist yet and PBS job $pbsJobId is not done; waiting...")
        Thread.sleep(5000)
        recursiveWait()
      case (false, true) =>
        logger.error(s"Output files do not exist and PBS job $pbsJobId is done")
         jobDescriptor.returnCode.clear().appendLine("69")
         69 // EX_UNAVAILABLE code to indicate PBS job is in Finished/Completed state, and output files not all created.
    }
    recursiveWait()
  }
  
  /**
   * Returns true if the PBS job is finished/completed; used to catch the case where job has been
   * killed or otherwise expired without creating rc file - e.g. if walltime was exceeded, or there
   * is a shell syntax error in user's command. 
   */
  private def jobDone(pbsJobId: Int): Boolean =
    if (installationIsPBSPro) {
      ( s"qstat -x -f ${pbsJobId}" #| Seq("grep", "job_state = F") ! ) == 0
    } else {
      ( s"qstat -f ${pbsJobId}" #| Seq("grep", "job_state = C") ! ) == 0
    }

  /**
   * This returns a zero-parameter function which, when called, will kill the
   * PBS job with id `pbsJobId`.  It also writes to the 'rc' file the value '143'
   */
  private def killPbsJob(jobDescriptor: BackendCallJobDescriptor, pbsJobId: Int) = () => {
    val qdelStdoutWriter = jobDescriptor.callRootPath.resolve("qdel.stdout").newBufferedWriter
    val qdelStderrWriter = jobDescriptor.callRootPath.resolve("qdel.stderr").newBufferedWriter
    val argv = Seq("qdel", pbsJobId.toString)
    val process = argv.run(ProcessLogger(qdelStdoutWriter writeWithNewline, qdelStderrWriter writeWithNewline))
    val returnCode: Int = process.exitValue()
    Vector(qdelStdoutWriter, qdelStderrWriter) foreach { _.flushAndClose() }
    jobDescriptor.returnCode.clear().appendLine("143")
    val logger = jobLogger(jobDescriptor)
    logger.debug(s"qdel $pbsJobId (returnCode=$returnCode)")
  }

  /**
   * Writes the script file containing the user's command from the WDL as well
   * as some extra shell code for monitoring jobs. Using a subprocess in the
   * script guarantees that rc file will be produced even if there are syntax
   * errors in the user's command.
   */
  private def writeScript(jobDescriptor: BackendCallJobDescriptor, instantiatedCommand: String) = {
    jobDescriptor.script.write(
      s"""#!/bin/sh
         |cd $$PBS_O_WORKDIR
         |sh -c 'set -e
         |${instantiatedCommand.replaceAll("'", "\"'\"")}
         |'
         |echo $$? > rc
         |""".stripMargin)
  }

  /**
   * Launches the qsub command, returns a tuple: (rc, Option(pbs_job_id))
   */
  private def launchQsub(jobDescriptor: BackendCallJobDescriptor): (Int, Option[Int]) = {
    val logger = jobLogger(jobDescriptor)
    val pbsJobName = s"crmwll_${jobDescriptor.workflowDescriptor.shortId}" take 15
    val queueSpec = jobDescriptor.callRuntimeAttributes.queue match {
      case Some(q) => List("-q", q)
      case None => List()
    }
    val argv = "qsub" :: queueSpec ++ Seq(
      "-N", pbsJobName,
      "-e", jobDescriptor.stderr.toAbsolutePath,
      "-o", jobDescriptor.stdout.toAbsolutePath,
      "-l", s"ncpus=${jobDescriptor.callRuntimeAttributes.cpu}",
      "-l", s"mem=${jobDescriptor.callRuntimeAttributes.memoryGB.toInt}gb",
      "-l", s"walltime=${jobDescriptor.callRuntimeAttributes.walltime}",
      jobDescriptor.script.toAbsolutePath
    ).map(_.toString)
    val backendCommandString = argv.map(s => "\""+s+"\"").mkString(" ")
    logger.info(s"backend command: $backendCommandString")

    val qsubStdoutFile = jobDescriptor.callRootPath.resolve("qsub.stdout")
    val qsubStderrFile = jobDescriptor.callRootPath.resolve("qsub.stderr")
    val qsubStdoutWriter = qsubStdoutFile.newBufferedWriter
    val qsubStderrWriter = qsubStderrFile.newBufferedWriter
    val executionDir = jobDescriptor.callRootPath.toAbsolutePath.toFile
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
  private def pollForPbsJobCompletionThenPostProcess(jobDescriptor: BackendCallJobDescriptor, pbsJobId: Int)(implicit ec: ExecutionContext): Future[(ExecutionResult, Int)] = {
    val logger = jobLogger(jobDescriptor)
    val abortFunction = killPbsJob(jobDescriptor, pbsJobId)
    val waitUntilCompleteFunction = waitUntilComplete(jobDescriptor, pbsJobId)
    jobDescriptor.abortRegistrationFunction.foreach(_.register(AbortFunction(abortFunction)))
    val jobReturnCode = waitUntilCompleteFunction
    val continueOnReturnCode = jobDescriptor.callRuntimeAttributes.continueOnReturnCode
    logger.info(s"PBS job completed (returnCode=$jobReturnCode)")
    val executionResult = (jobReturnCode, jobDescriptor.stderr.toFile.length) match {
      case (r, _) if r == 143 => AbortedExecution.future // Special case to check for SIGTERM exit code - implying abort
      case (r, _) if r == 69 => // Special case to check for PBS job not appearing in qstat
        val message = s"PBS job $pbsJobId reached defunct state without all output files being created"
        logger.error(message)
        NonRetryableExecution(new Exception(message), Option(r)).future
      case (r, _) if !continueOnReturnCode.continueFor(r) =>
        val message = s"PBS job failed because of return code: $r"
        logger.error(message)
        NonRetryableExecution(new Exception(message), Option(r)).future
      case (_, stderrLength) if stderrLength > 0 && jobDescriptor.callRuntimeAttributes.failOnStderr =>
        val message = s"PBS job failed because there were $stderrLength bytes on standard error"
        logger.error(message)
        NonRetryableExecution(new Exception(message), Option(0)).future
      case (r, _) =>
        postProcess(jobDescriptor) match {
          case Success(callOutputs) => jobDescriptor.hash map { h => SuccessfulBackendCallExecution(callOutputs, Seq.empty, r, h) }
          case Failure(e) => NonRetryableExecution(e).future
        }
    }
    executionResult map { (_,  jobReturnCode) }
  }

  override def callRootPathWithBaseRoot(descriptor: BackendCallJobDescriptor, baseRoot: String): Path = {
    val path = super.callRootPathWithBaseRoot(descriptor, baseRoot)
    if (!path.toFile.exists()) path.toFile.mkdirs()
    path
  }

  override def executionInfoKeys: List[String] = List(InfoKeys.JobNumber)

  override def callEngineFunctions(descriptor: BackendCallJobDescriptor): CallEngineFunctions = {
    new PbsCallEngineFunctions(descriptor.workflowDescriptor.fileSystems, buildCallContext(descriptor))
  }

  override def fileSystems(options: WorkflowOptions): List[FileSystem] = {
    val gcsStorage = StorageFactory.userAuthenticated(options) orElse StorageFactory.cromwellAuthenticated
    val gcs = gcsStorage map GcsFileSystemProvider.apply map { _.getFileSystem } toOption

    List(gcs, Option(defaultFileSystem)).flatten
  }

  def instantiateCommand(jobDescriptor: BackendCallJobDescriptor): Try[String] = {
    val backendInputs = adjustInputPaths(jobDescriptor)
    jobDescriptor.call.instantiateCommandLine(backendInputs, jobDescriptor.callEngineFunctions)
  }

  override def poll(jobDescriptor: BackendCallJobDescriptor, previous: ExecutionHandle)(implicit ec: ExecutionContext) = Future.successful(previous)

  override def resume(descriptor: BackendCallJobDescriptor, jobKey: BackendJobKey)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    Future.failed(new Throwable("resume invoked on non-resumable PBS backend"))
  }
}
