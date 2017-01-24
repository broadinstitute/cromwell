package cromwell.backend.sfs

import java.nio.file.{FileAlreadyExistsException, Path}

import better.files._
import cromwell.backend._
import cromwell.backend.async.{ExecutionHandle, FailedNonRetryableExecutionHandle, PendingExecutionHandle}
import cromwell.backend.io.JobPathsWithDocker
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncJob}
import cromwell.backend.validation._
import cromwell.core.path.DefaultPathBuilder
import cromwell.core.path.FileImplicits._
import cromwell.core.path.PathFactory._
import cromwell.core.retry.SimpleExponentialBackoff
import wdl4s.values.WdlFile

import scala.concurrent.duration._

case class SharedFileSystemRunStatus(returnCodeFileExists: Boolean) {
  override def toString: String = if (returnCodeFileExists) "Done" else "WaitingForReturnCodeFile"
}

object SharedFileSystemAsyncJobExecutionActor {
  val JobIdKey = "sfs_job_id"
}

/**
  * Runs a job on a shared backend, with the ability to (abstractly) submit asynchronously, then poll, kill, etc.
  *
  * Abstract workhorse of the shared file system.
  *
  * The trait requires that there exist:
  * - Some unix process to submit jobs asynchronously.
  * - When the job runs, outputs will be written to the same filesystem cromwell is executing.
  *
  * As the job runs, this backend will poll for an `rc` file. The `rc` file should be written after the command
  * completes, or will be written by this trait itself during an abort.
  *
  * In practice instead of extending this trait, most systems requiring a backend can likely just configure a backend in
  * the application.conf using a cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory
  *
  *
  * NOTE: Although some methods return futures due to the (current) contract in BJEA/ABJEA, this actor only executes
  * during the receive, and does not launch new runnables/futures from inside "receive"... except--
  *
  * The __one__ exception is that the when `poll` is processing a successful return code. Currently processReturnCode
  * is calling into a stub for generating fake hashes. This functionality is TBD, but it is likely that we __should__
  * begin teardown while we return a future with the results, assuming we're still using futures instead of akka-ish
  * messages.
  */
trait SharedFileSystemAsyncJobExecutionActor
  extends BackendJobLifecycleActor with StandardAsyncExecutionActor with SharedFileSystemJobCachingActorHelper {

  override type StandardAsyncRunInfo = Any

  override type StandardAsyncRunStatus = SharedFileSystemRunStatus

  override lazy val pollBackOff = SimpleExponentialBackoff(1.second, 5.minutes, 1.1)

  override lazy val executeOrRecoverBackOff = SimpleExponentialBackoff(3.seconds, 30.seconds, 1.1)

  /**
    * Returns the command for running the job. The returned command may or may not run the job asynchronously in the
    * background. If the command does not run the script asynchronously in the background or on some job scheduler, the
    * trait `BackgroundAsyncJobExecutionActor` should be mixed in to run these processArgs inside a bash script in the
    * background.
    *
    * @return The command to run a script.
    */
  def processArgs: SharedFileSystemCommand

  /**
    * Retrieves the job id after the command has been submited for asynchronous running.
    *
    * @param exitValue The exit value of the submit.
    * @param stdout    The stdout of the submit.
    * @param stderr    The stderr of the submit.
    * @return The job id wrapped in a SharedFileSystemJob.
    */
  def getJob(exitValue: Int, stdout: Path, stderr: Path): StandardAsyncJob

  /**
    * Returns the command for checking if a job is alive, returing non-zero if the job cannot be found or has errored.
    *
    * @param job The job to check.
    * @return The command for checking if a job is alive.
    */
  def checkAliveArgs(job: StandardAsyncJob): SharedFileSystemCommand

  /**
    * Returns the command for killing a job.
    *
    * @param job The job to kill.
    * @return The command for killing a job.
    */
  def killArgs(job: StandardAsyncJob): SharedFileSystemCommand

  lazy val jobPathsWithDocker: JobPathsWithDocker = jobPaths.asInstanceOf[JobPathsWithDocker]

  def jobName: String = s"cromwell_${jobDescriptor.workflowDescriptor.id.shortString}_${jobDescriptor.call.unqualifiedName}"

  lazy val isDockerRun: Boolean = RuntimeAttributesValidation.extractOption(
    DockerValidation.instance, validatedRuntimeAttributes).isDefined

  /**
    * Localizes the file, run outside of docker.
    */
  override def preProcessWdlFile(wdlFile: WdlFile): WdlFile = {
    sharedFileSystem.localizeWdlFile(jobPathsWithDocker.callInputsRoot, isDockerRun)(wdlFile)
  }

  /**
    * Returns the paths to the file, inside of docker.
    */
  override def mapCommandLineWdlFile(wdlFile: WdlFile): WdlFile = {
    val cleanPath = DefaultPathBuilder.build(wdlFile.valueString).get
    WdlFile(if (isDockerRun) jobPathsWithDocker.toDockerPath(cleanPath).toString else cleanPath.toString)
  }

  override lazy val commandDirectory: Path = {
    if (isDockerRun) jobPathsWithDocker.callExecutionDockerRoot else jobPaths.callExecutionRoot
  }

  override def execute(): ExecutionHandle = {
    File(jobPaths.callExecutionRoot).createPermissionedDirectories()
    writeScriptContents()
    val runner = makeProcessRunner()
    val exitValue = runner.run()
    if (exitValue != 0) {
      FailedNonRetryableExecutionHandle(new RuntimeException("Unable to start job. " +
        s"Check the stderr file for possible errors: ${runner.stderrPath}"))
    } else {
      val runningJob = getJob(exitValue, runner.stdoutPath, runner.stderrPath)
      PendingExecutionHandle(jobDescriptor, runningJob, None, None)
    }
  }

  def writeScriptContents(): Unit = {
    File(jobPaths.script).write(commandScriptContents)
    ()
  }

  /**
    * Creates a script to submit the script for asynchronous processing. The default implementation assumes the
    * processArgs already runs the script asynchronously. If not, mix in the `BackgroundAsyncJobExecutionActor` that
    * will run the command in the background, and return a PID for the backgrounded process.
    *
    * @return A process runner that will relatively quickly submit the script asynchronously.
    */
  def makeProcessRunner(): ProcessRunner = {
    val stdout = pathPlusSuffix(jobPaths.stdout, "submit")
    val stderr = pathPlusSuffix(jobPaths.stderr, "submit")
    new ProcessRunner(processArgs.argv, stdout.path, stderr.path)
  }

  override def recover(job: StandardAsyncJob): ExecutionHandle = {
    // To avoid race conditions, check for the rc file after checking if the job is alive.
    if (isAlive(job) || File(jobPaths.returnCode).exists) {
      // If we're done, we'll get to the rc during the next poll.
      // Or if we're still running, return pending also.
      jobLogger.info(s"Recovering using job id: ${job.jobId}")
      PendingExecutionHandle(jobDescriptor, job, None, None)
    } else {
      // Could start executeScript(), but for now fail because we shouldn't be in this state.
      FailedNonRetryableExecutionHandle(new RuntimeException(
        s"Unable to determine that ${job.jobId} is alive, and ${jobPaths.returnCode} does not exist."), None)
    }
  }

  def isAlive(job: StandardAsyncJob): Boolean = {
    val argv = checkAliveArgs(job).argv
    val stdout = pathPlusSuffix(jobPaths.stdout, "check")
    val stderr = pathPlusSuffix(jobPaths.stderr, "check")
    val checkAlive = new ProcessRunner(argv, stdout.path, stderr.path)
    checkAlive.run() == 0
  }

  override def tryAbort(job: StandardAsyncJob): Unit = {
    val returnCodeTmp = pathPlusSuffix(jobPaths.returnCode, "kill")
    returnCodeTmp.write(s"$SIGTERM\n")
    try {
      returnCodeTmp.moveTo(jobPaths.returnCode)
    } catch {
      case _: FileAlreadyExistsException =>
        // If the process has already completed, there will be an existing rc file.
        returnCodeTmp.delete(true)
    }
    val argv = killArgs(job).argv
    val stdout = pathPlusSuffix(jobPaths.stdout, "kill")
    val stderr = pathPlusSuffix(jobPaths.stderr, "kill")
    val killer = new ProcessRunner(argv, stdout.path, stderr.path)
    killer.run()
    ()
  }

  override def pollStatus(handle: StandardAsyncPendingExecutionHandle): SharedFileSystemRunStatus = {
    SharedFileSystemRunStatus(File(jobPaths.returnCode).exists)
  }

  override def isTerminal(runStatus: StandardAsyncRunStatus): Boolean = {
    runStatus.returnCodeFileExists
  }

  override def mapOutputWdlFile(wdlFile: WdlFile): WdlFile = {
    sharedFileSystem.mapJobWdlFile(jobPaths)(wdlFile)
  }
}
