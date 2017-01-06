package cromwell.backend.sfs

import java.nio.file.{FileAlreadyExistsException, Path}

import akka.actor.{Actor, ActorLogging}
import better.files._
import cromwell.backend._
import cromwell.backend.async.{AsyncBackendJobExecutionActor, ExecutionHandle, FailedNonRetryableExecutionHandle, PendingExecutionHandle, SuccessfulExecutionHandle}
import cromwell.backend.io.WorkflowPathsBackendInitializationData
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncJob}
import cromwell.backend.validation._
import cromwell.backend.wdl.OutputEvaluator
import cromwell.core.WorkflowId
import cromwell.core.path.PathFactory._
import cromwell.core.path.{DefaultPathBuilder, PathBuilder}
import cromwell.core.retry.SimpleExponentialBackoff
import wdl4s.values.{WdlArray, WdlFile, WdlGlobFile, WdlMap, WdlValue}
import wdl4s.{EvaluatedTaskInputs, TaskCall}

import scala.concurrent.duration._
import scala.util.{Failure, Success, Try}

object SharedFileSystemJob {
  val JobIdKey = "sfs_job_id"
}

case class SharedFileSystemRunStatus(returnCodeFileExists: Boolean)

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
  extends Actor with ActorLogging with BackendJobLifecycleActor with AsyncBackendJobExecutionActor
    with StandardAsyncExecutionActor with SharedFileSystemJobCachingActorHelper {

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

  def toUnixPath(docker: Boolean)(path: WdlValue): WdlValue = {
    path match {
      case _: WdlFile =>
        val cleanPath = DefaultPathBuilder.build(path.valueString).get
        WdlFile(if (docker) jobPaths.toDockerPath(cleanPath).toString else cleanPath.toString)
      case array: WdlArray => WdlArray(array.wdlType, array.value map toUnixPath(docker))
      case map: WdlMap => WdlMap(map.wdlType, map.value mapValues toUnixPath(docker))
      case wdlValue => wdlValue
    }
  }

  def jobName: String = s"cromwell_${jobDescriptor.workflowDescriptor.id.shortString}_${jobDescriptor.call.unqualifiedName}"

  lazy val workflowDescriptor: BackendWorkflowDescriptor = jobDescriptor.workflowDescriptor
  lazy val call: TaskCall = jobDescriptor.key.call
  lazy val pathBuilders: List[PathBuilder] = WorkflowPathsBackendInitializationData.pathBuilders(backendInitializationDataOption)
  private[sfs] lazy val backendEngineFunctions = SharedFileSystemExpressionFunctions(jobPaths, pathBuilders)
  override lazy val workflowId: WorkflowId = jobDescriptor.workflowDescriptor.id
  override lazy val jobTag: String = jobDescriptor.key.tag

  lazy val isDockerRun: Boolean = RuntimeAttributesValidation.extractOption(
    DockerValidation.instance, validatedRuntimeAttributes).isDefined

  override lazy val commandLineFunctions: SharedFileSystemExpressionFunctions = backendEngineFunctions

  override lazy val commandLinePreProcessor: (EvaluatedTaskInputs) => Try[EvaluatedTaskInputs] =
    sharedFileSystem.localizeInputs(jobPaths.callInputsRoot, isDockerRun)

  override lazy val commandLineValueMapper: (WdlValue) => WdlValue = toUnixPath(isDockerRun)

  override lazy val startMetadataKeyValues: Map[String, Any] = {
    super[SharedFileSystemJobCachingActorHelper].startMetadataKeyValues
  }

  override def execute(): ExecutionHandle = {
    val script = instantiatedCommand
    jobLogger.info(s"`$script`")
    File(jobPaths.callExecutionRoot).createDirectories()
    val cwd = if (isDockerRun) jobPaths.callExecutionDockerRoot else jobPaths.callExecutionRoot
    writeScript(script, cwd, backendEngineFunctions.findGlobOutputs(call, jobDescriptor))
    jobLogger.info(s"command: $processArgs")
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

  /**
    * Writes the script file containing the user's command from the WDL as well
    * as some extra shell code for monitoring jobs
    */
  private def writeScript(instantiatedCommand: String, cwd: Path, globFiles: Set[WdlGlobFile]) = {
    val rcPath = if (isDockerRun) jobPaths.toDockerPath(jobPaths.returnCode) else jobPaths.returnCode
    val rcTmpPath = pathPlusSuffix(rcPath, "tmp").path

    def globManipulation(globFile: WdlGlobFile) = {

      // TODO: Move glob list and directory generation into trait GlobFunctions? There is already a globPath using callContext
      val globDir = backendEngineFunctions.globName(globFile.value)
      val globDirectory = File(cwd)./(globDir)
      val globList = File(cwd)./(s"$globDir.list")

      s"""|mkdir $globDirectory
          |( ln -L ${globFile.value} $globDirectory 2> /dev/null ) || ( ln ${globFile.value} $globDirectory )
          |ls -1 $globDirectory > $globList
          |""".stripMargin
    }

    val globManipulations = globFiles.map(globManipulation).mkString("\n")

    val scriptBody =
      s"""|#!/bin/sh
          |(
          |cd $cwd
          |INSTANTIATED_COMMAND
          |)
          |echo $$? > $rcTmpPath
          |(
          |cd $cwd
          |$globManipulations
          |)
          |mv $rcTmpPath $rcPath
          |""".stripMargin.replace("INSTANTIATED_COMMAND", instantiatedCommand)

    File(jobPaths.script).write(scriptBody)
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

  override def remoteStdErrPath: Path = jobPaths.stderr

  override def remoteReturnCodePath: Path = jobPaths.returnCode

  override def continueOnReturnCode: ContinueOnReturnCode = RuntimeAttributesValidation.extract(
    ContinueOnReturnCodeValidation.instance, validatedRuntimeAttributes)

  override def failOnStdErr: Boolean = RuntimeAttributesValidation.extract(
    FailOnStderrValidation.instance, validatedRuntimeAttributes)

  override def pollStatus(handle: StandardAsyncPendingExecutionHandle): SharedFileSystemRunStatus = {
    SharedFileSystemRunStatus(File(jobPaths.returnCode).exists)
  }

  override def isTerminal(runStatus: StandardAsyncRunStatus): Boolean = {
    runStatus.returnCodeFileExists
  }

  override def handleExecutionSuccess(runStatus: StandardAsyncRunStatus, handle: StandardAsyncPendingExecutionHandle,
                                      returnCode: Int): ExecutionHandle = {
    val outputsTry =
      OutputEvaluator.evaluateOutputs(jobDescriptor, backendEngineFunctions, sharedFileSystem.outputMapper(jobPaths))
    outputsTry match {
      case Success(outputs) => SuccessfulExecutionHandle(outputs, returnCode, jobPaths.detritusPaths, Seq.empty)
      case Failure(throwable) => FailedNonRetryableExecutionHandle(throwable, Option(returnCode))
    }
  }

}
