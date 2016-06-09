package cromwell.engine.backend.lsf

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
import cromwell.engine.backend.lsf.LsfBackend.InfoKeys
import cromwell.engine.db.DataAccess._
import cromwell.engine.workflow.BackendCallKey
import cromwell.logging.WorkflowLogger
import cromwell.util.FileUtil._

import com.typesafe.config.ConfigFactory
import scala.annotation.tailrec
import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.sys.process._
import scala.util.{Failure, Success, Try}
import scala.collection.JavaConversions._

object LsfBackend {
  object InfoKeys {
    val JobNumber = "LSF_JOB_NUMBER"
  }

  implicit class LsfEnhancedJobDescriptor(val jobDescriptor: BackendCallJobDescriptor) extends AnyVal {
    def workflowRootPath = jobDescriptor.workflowDescriptor.workflowRootPath
    def stdout = jobDescriptor.callRootPath.resolve("stdout")
    def stderr = jobDescriptor.callRootPath.resolve("stderr")
    def script = jobDescriptor.callRootPath.resolve("script.sh")
    def returnCode = jobDescriptor.callRootPath.resolve("rc")
  }
}

case class LsfBackend(actorSystem: ActorSystem) extends Backend with SharedFileSystem {
  def returnCode(jobDescriptor: BackendCallJobDescriptor) = jobDescriptor.returnCode

  import LocalBackend.WriteWithNewline

  override def backendType = BackendType.LSF

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
          case (bsubReturnCode, _) if bsubReturnCode != 0 => NonRetryableExecution(
            new Throwable(s"Error: bsub exited with return code: $bsubReturnCode")).future
          case (_, None) => NonRetryableExecution(new Throwable(s"Could not parse Job ID from bsub output")).future
          case (_, Some(lsfJobId)) =>
            val updateDatabaseWithRunningInfo = updateLsfJobTable(jobDescriptor, "Running", None, Option(lsfJobId))
            pollForLsfJobCompletionThenPostProcess(jobDescriptor, lsfJobId) map { case (executionResult, jobRc) =>
              // Only send the completion update once the 'Running' update has completed (regardless of the success of that update)
              updateDatabaseWithRunningInfo onComplete {
                case _ =>
                  val completionStatus = statusString(executionResult)
                  val updateDatabaseWithCompletionInfo = updateLsfJobTable(jobDescriptor, completionStatus, Option(jobRc), Option(lsfJobId))
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

  private def updateLsfJobTable(jobDescriptor: BackendCallJobDescriptor, status: String, rc: Option[Int], lsfJobId: Option[Int])
                               (implicit ec: ExecutionContext): Future[Unit] = {
    globalDataAccess.updateExecutionInfo(jobDescriptor.workflowDescriptor.id, BackendCallKey(jobDescriptor.call, jobDescriptor.key.index, jobDescriptor.key.attempt), InfoKeys.JobNumber, Option(lsfJobId.toString))
  }

  /** TODO restart isn't currently implemented for LSF, there is probably work that needs to be done here much like
    * JES restart, which perhaps could be factored out into a common "remote executor" trait.
    */
  override def prepareForRestart(restartableWorkflow: WorkflowDescriptor)(implicit ec: ExecutionContext) = Future.successful(())

  /**
   * Returns the RC of this job when it finishes.  Sleeps and polls
   * until the 'rc' file is generated
   */
  private def waitUntilComplete(jobDescriptor: BackendCallJobDescriptor): Int = {
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
   * LSF job with id `lsfJobId`.  It also writes to the 'rc' file the value '143'
   */
  private def killLsfJob(jobDescriptor: BackendCallJobDescriptor, lsfJobId: Int) = () => {
    val bkillStdoutWriter = jobDescriptor.callRootPath.resolve("bkill.stdout").newBufferedWriter
    val bkillStderrWriter = jobDescriptor.callRootPath.resolve("bkill.stderr").newBufferedWriter
    val argv = Seq("bkill", lsfJobId.toString)
    val process = argv.run(ProcessLogger(bkillStdoutWriter writeWithNewline, bkillStderrWriter writeWithNewline))
    val returnCode: Int = process.exitValue()
    Vector(bkillStdoutWriter, bkillStderrWriter) foreach { _.flushAndClose() }
    jobDescriptor.returnCode.clear().appendLine("143")
    val logger = jobLogger(jobDescriptor)
    logger.debug(s"bkill $lsfJobId (returnCode=$returnCode)")
  }

  /**
   * Writes the script file containing the user's command from the WDL as well
   * as some extra shell code for monitoring jobs
   */
  private def writeScript(jobDescriptor: BackendCallJobDescriptor, instantiatedCommand: String) = {
    jobDescriptor.script.write(
      s"""#!/bin/sh
         |$instantiatedCommand
         |echo $$? > rc
         |""".stripMargin)
  }

  /**
   * Launches the bsub command, returns a tuple: (rc, Option(lsf_job_id))
   */
  private def launchQsub(jobDescriptor: BackendCallJobDescriptor): (Int, Option[Int]) = {
    val logger = jobLogger(jobDescriptor)
    val backendConf = ConfigFactory.load.getConfig("backend")
    val lsfConf = backendConf.getConfig("lsf")
 
    val lsfOption = lsfConf.root().unwrapped()
    if(lsfOption.get("-J") == null) {
       lsfOption.put("-J", s"cromwell_${jobDescriptor.workflowDescriptor.shortId}_${jobDescriptor.call.unqualifiedName}")
    }
    if(lsfOption.get("-cwd") == null) {
       lsfOption.put("-cwd", jobDescriptor.callRootPath.toAbsolutePath)
    }
    if(lsfOption.get("-o") == null) {
       lsfOption.put("-o", jobDescriptor.stdout.getFileName)
    }
    if(lsfOption.get("-e") == null) {
       lsfOption.put("-e", jobDescriptor.stderr.getFileName)
    }

    val argv = (lsfOption.foldLeft(Seq("bsub"))((command, kv) =>  command ++ Seq(kv._1, kv._2.toString)) ++  Seq("/bin/sh", jobDescriptor.script.toAbsolutePath)).map(_.toString)
    val backendCommandString = argv.map(s => "\""+s+"\"").mkString(" ")
    logger.info(s"backend command: $backendCommandString")

    val bsubStdoutFile = jobDescriptor.callRootPath.resolve("bsub.stdout")
    val bsubStderrFile = jobDescriptor.callRootPath.resolve("bsub.stderr")
    val bsubStdoutWriter = bsubStdoutFile.newBufferedWriter
    val bsubStderrWriter = bsubStderrFile.newBufferedWriter
    val process = argv.run(ProcessLogger(bsubStdoutWriter writeWithNewline, bsubStderrWriter writeWithNewline))
    val returnCode: Int = process.exitValue()
    Vector(bsubStdoutWriter, bsubStderrWriter) foreach { _.flushAndClose() }

    // The bsub makes it so stdout only has the job ID, if it was successfully launched
    val pattern = "Job <(\\d+)>".r;
    val jobId = Try(pattern.findFirstIn(bsubStdoutFile.contentAsString.stripLineEnd)) match {
      case Success(line) => line.getOrElse("0") match { case pattern(id) => Some(id.toInt) case _ => Some(0) }
      case Failure(ex) =>
        logger.error(s"Could not find LSF job ID from bsub stdout file.\n\n" +
          s"Check the bsub stderr file for possible errors: ${bsubStderrFile.toAbsolutePath}")
        None
    }
    logger.info(s"bsub returnCode=$returnCode, job ID=${jobId.getOrElse("NONE")}")
    (returnCode, jobId)
  }

  /**
   * This waits for a given LSF job to finish.  When finished, it post-processes the job
   * and returns the outputs for the jobDescriptor
   */
  private def pollForLsfJobCompletionThenPostProcess(jobDescriptor: BackendCallJobDescriptor, lsfJobId: Int)(implicit ec: ExecutionContext): Future[(ExecutionResult, Int)] = {
    val logger = jobLogger(jobDescriptor)
    val abortFunction = killLsfJob(jobDescriptor, lsfJobId)
    val waitUntilCompleteFunction = waitUntilComplete(jobDescriptor)
    jobDescriptor.abortRegistrationFunction.foreach(_.register(AbortFunction(abortFunction)))
    val jobReturnCode = waitUntilCompleteFunction
    val continueOnReturnCode = jobDescriptor.callRuntimeAttributes.continueOnReturnCode
    logger.info(s"LSF job completed (returnCode=$jobReturnCode)")
    val executionResult = (jobReturnCode, jobDescriptor.stderr.toFile.length) match {
      case (r, _) if r == 143 => AbortedExecution.future // Special case to check for SIGTERM exit code - implying abort
      case (r, _) if !continueOnReturnCode.continueFor(r) =>
        val message = s"LSF job failed because of return code: $r"
        logger.error(message)
        NonRetryableExecution(new Exception(message), Option(r)).future
      case (_, stderrLength) if stderrLength > 0 && jobDescriptor.callRuntimeAttributes.failOnStderr =>
        val message = s"LSF job failed because there were $stderrLength bytes on standard error"
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
    new LsfCallEngineFunctions(descriptor.workflowDescriptor.fileSystems, buildCallContext(descriptor))
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
    Future.failed(new Throwable("resume invoked on non-resumable LSF backend"))
  }
}
