package cromwell.backend.sfs

import java.nio.file.FileAlreadyExistsException
import java.time.Instant

import cromwell.backend._
import cromwell.backend.async.{ExecutionHandle, FailedNonRetryableExecutionHandle, PendingExecutionHandle}
import cromwell.backend.io.JobPathsWithDocker
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncJob}
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.core.retry.SimpleExponentialBackoff
import wom.values.WomFile

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.util.{Failure, Success, Try}

sealed trait SharedFileSystemRunState {
  def status: String
  def terminal: Boolean

  override def toString: String = status
}

case class SharedFileSystemJobRunning(validUntil: Option[Instant]) extends SharedFileSystemRunState {
  override def terminal: Boolean = false
  override def status = "Running"

  // Whether this running state is stale (ie has the 'validUntil' time passed?)
  def stale: Boolean = validUntil.exists(t => t.isBefore(Instant.now))
}

case class SharedFileSystemJobWaitingForReturnCode(waitUntil: Option[Instant]) extends SharedFileSystemRunState {
  override def terminal: Boolean = false
  override def status = "WaitingForReturnCode"

  // Whether or not to give up waiting for the RC to appear (ie has the 'waitUntil' time passed?)
  def giveUpWaiting: Boolean = waitUntil.exists(_.isBefore(Instant.now))
}

case object SharedFileSystemJobDone extends SharedFileSystemRunState {
  override def terminal: Boolean = true
  override def status = "Done"
}
case object SharedFileSystemJobFailed extends SharedFileSystemRunState {
  override def terminal: Boolean = true
  override def status = "Failed"
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

  override type StandardAsyncRunState = SharedFileSystemRunState

  /** True if the status contained in `thiz` is equivalent to `that`, delta any other data that might be carried around
    * in the state type.
    */
  override def statusEquivalentTo(thiz: StandardAsyncRunState)(that: StandardAsyncRunState): Boolean = thiz.status == that.status

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

  def jobName: String = s"cromwell_${jobDescriptor.workflowDescriptor.id.shortString}_${jobDescriptor.taskCall.localName}"

  /**
    * Localizes the file, run outside of docker.
    */
  override def preProcessWomFile(womFile: WomFile): WomFile = {
    sharedFileSystem.localizeWomFile(jobPathsWithDocker.callInputsRoot, isDockerRun)(womFile)
  }

  /**
    * Returns the paths to the file, inside of docker.
    */
  override def mapCommandLineWomFile(womFile: WomFile): WomFile = {
    womFile mapFile { path =>
      val cleanPath = DefaultPathBuilder.build(path).get
      if (isDockerRun) jobPathsWithDocker.toDockerPath(cleanPath).pathAsString else cleanPath.pathAsString
    }
  }

  override lazy val commandDirectory: Path = {
    if (isDockerRun) jobPathsWithDocker.callExecutionDockerRoot else jobPaths.callExecutionRoot
  }

  override def execute(): ExecutionHandle = {
    if (isDockerRun) jobPaths.callExecutionRoot.createPermissionedDirectories()
    else jobPaths.callExecutionRoot.createDirectories()
    writeScriptContents().fold(
      identity,
      { _ =>
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
    )
  }

  def writeScriptContents(): Either[ExecutionHandle, Unit] =
    commandScriptContents.fold(
      errors => Left(FailedNonRetryableExecutionHandle(new RuntimeException("Unable to start job due to: " + errors.toList.mkString(", ")))),
      {script => jobPaths.script.write(script); Right(())} )

  lazy val standardPaths = jobPaths.standardPaths
  /**
    * Creates a script to submit the script for asynchronous processing. The default implementation assumes the
    * processArgs already runs the script asynchronously. If not, mix in the `BackgroundAsyncJobExecutionActor` that
    * will run the command in the background, and return a PID for the backgrounded process.
    *
    * @return A process runner that will relatively quickly submit the script asynchronously.
    */
  def makeProcessRunner(): ProcessRunner = {
    val stdout = standardPaths.output.plusExt("submit")
    val stderr = standardPaths.error.plusExt("submit")
    new ProcessRunner(processArgs.argv, stdout, stderr)
  }

  override def recover(job: StandardAsyncJob): ExecutionHandle = reconnectToExistingJob(job)

  override def reconnectAsync(job: StandardAsyncJob): Future[ExecutionHandle] = {
    Future.successful(reconnectToExistingJob(job))
  }

  override def reconnectToAbortAsync(job: StandardAsyncJob): Future[ExecutionHandle] = {
    Future.successful(reconnectToExistingJob(job, forceAbort = true))
  }

  private def reconnectToExistingJob(job: StandardAsyncJob, forceAbort: Boolean = false): ExecutionHandle = {
    // To avoid race conditions, check for the rc file after checking if the job is alive.
    isAlive(job) match {
      case Success(true) =>
        // If the job is not done and forceAbort is true, try to abort it
        if (!jobPaths.returnCode.exists && forceAbort) {
          jobLogger.info(s"Recovering and aborting using job id: ${job.jobId}")
          tryAbort(job)
        }
        PendingExecutionHandle(jobDescriptor, job, None, None)
      case Success(false) =>
        if (jobPaths.returnCode.exists) {
          PendingExecutionHandle(jobDescriptor, job, None, None)
        } else {
          // Could start executeScript(), but for now fail because we shouldn't be in this state.
          FailedNonRetryableExecutionHandle(new RuntimeException(
            s"Unable to determine that ${job.jobId} is alive, and ${jobPaths.returnCode} does not exist."), None)
        }
      case Failure(f) => FailedNonRetryableExecutionHandle(f, None)
    }
  }

  def isAlive(job: StandardAsyncJob): Try[Boolean] = Try {
    val argv = checkAliveArgs(job).argv
    val stdout = standardPaths.output.plusExt("check")
    val stderr = standardPaths.error.plusExt("check")
    val checkAlive = new ProcessRunner(argv, stdout, stderr)
    checkAlive.run() == 0
  }

  override def requestsAbortAndDiesImmediately: Boolean = false

  override def tryAbort(job: StandardAsyncJob): Unit = {
    val returnCodeTmp = jobPaths.returnCode.plusExt("kill")
    returnCodeTmp.write(s"$SIGTERM\n")
    try {
      returnCodeTmp.moveTo(jobPaths.returnCode)
    } catch {
      case _: FileAlreadyExistsException =>
        // If the process has already completed, there will be an existing rc file.
        returnCodeTmp.delete(true)
    }
    val stderrTmp = standardPaths.error.plusExt("kill")
    stderrTmp.touch()
    try {
      stderrTmp.moveTo(standardPaths.error)
    } catch {
      case _: FileAlreadyExistsException =>
        // If the process has already started, there will be an existing stderr file.
        stderrTmp.delete(true)
    }
    val argv = killArgs(job).argv
    val stdout = standardPaths.output.plusExt("kill")
    val stderr = standardPaths.error.plusExt("kill")
    val killer = new ProcessRunner(argv, stdout, stderr)
    killer.run()
    ()
  }

  override def pollStatus(handle: StandardAsyncPendingExecutionHandle): SharedFileSystemRunState = {
    if (jobPaths.returnCode.exists) SharedFileSystemJobDone
    else SharedFileSystemJobWaitingForReturnCode(waitUntil = None)
  }

  override def isTerminal(runStatus: StandardAsyncRunState): Boolean = {
    runStatus.terminal
  }

  override def mapOutputWomFile(womFile: WomFile): WomFile = {
    sharedFileSystem.mapJobWomFile(jobPaths)(womFile)
  }
}
