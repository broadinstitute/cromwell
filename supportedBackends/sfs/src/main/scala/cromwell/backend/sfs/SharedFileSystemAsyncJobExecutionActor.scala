package cromwell.backend.sfs

import java.nio.file.{FileAlreadyExistsException, Path}

import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.event.LoggingReceive
import better.files._
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.async.AsyncBackendJobExecutionActor._
import cromwell.backend.async.{AbortedExecutionHandle, AsyncBackendJobExecutionActor, ExecutionHandle, FailedNonRetryableExecutionHandle, NonRetryableExecution, SuccessfulExecutionHandle}
import cromwell.backend.io.WorkflowPathsBackendInitializationData
import cromwell.backend.sfs.SharedFileSystem._
import cromwell.backend.validation._
import cromwell.backend.wdl.{OutputEvaluator, Command}
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptor}
import cromwell.core.JobOutputs
import cromwell.core.logging.JobLogging
import cromwell.core.path.DefaultPathBuilder
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.services.keyvalue.KeyValueServiceActor._
import wdl4s.values.{WdlArray, WdlFile, WdlMap, WdlValue}

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.util.{Failure, Success, Try}

object SharedFileSystemJob {
  val JobIdKey = "sfs_job_id"
}

/**
  * A generic job that runs and tracks some string identifier for the job.
  */
case class SharedFileSystemJob(jobId: String) extends JobId

case class SharedFileSystemAsyncJobExecutionActorParams
(
  serviceRegistryActor: ActorRef,
  jobDescriptor: BackendJobDescriptor,
  configurationDescriptor: BackendConfigurationDescriptor,
  completionPromise: Promise[BackendJobExecutionResponse],
  backendInitializationDataOption: Option[BackendInitializationData]
)

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
  extends Actor with ActorLogging with AsyncBackendJobExecutionActor with SharedFileSystemJobCachingActorHelper
    with JobLogging {

  case class SharedFileSystemPendingExecutionHandle(jobDescriptor: BackendJobDescriptor,
                                                    run: SharedFileSystemJob) extends ExecutionHandle {
    override val isDone = false
    override val result = NonRetryableExecution(new IllegalStateException(
      "SharedFileSystemPendingExecutionHandle cannot yield a result"))
  }

  context.become(sharedReceive(None) orElse super.receive)

  val SIGTERM = 143
  val SIGINT = 130

  override lazy val pollBackOff = SimpleExponentialBackoff(1.second, 5.minutes, 1.1)

  override lazy val executeOrRecoverBackOff = SimpleExponentialBackoff(3.seconds, 30.seconds, 1.1)

  override protected implicit def ec = context.dispatcher

  val params: SharedFileSystemAsyncJobExecutionActorParams

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
  def getJob(exitValue: Int, stdout: Path, stderr: Path): SharedFileSystemJob

  /**
    * Returns the command for checking if a job is alive, returing non-zero if the job cannot be found or has errored.
    *
    * @param job The job to check.
    * @return The command for checking if a job is alive.
    */
  def checkAliveArgs(job: SharedFileSystemJob): SharedFileSystemCommand

  /**
    * Returns the command for killing a job.
    *
    * @param job The job to kill.
    * @return The command for killing a job.
    */
  def killArgs(job: SharedFileSystemJob): SharedFileSystemCommand

  override lazy val jobDescriptor = params.jobDescriptor

  override lazy val completionPromise = params.completionPromise

  override lazy val serviceRegistryActor = params.serviceRegistryActor

  override lazy val configurationDescriptor = params.configurationDescriptor

  override lazy val backendInitializationDataOption = params.backendInitializationDataOption

  def toUnixPath(docker: Boolean)(path: WdlValue): WdlValue = {
    path match {
      case file: WdlFile =>
        val cleanPath = DefaultPathBuilder.build(path.valueString).get
        WdlFile(if (docker) jobPaths.toDockerPath(cleanPath).toString else cleanPath.toString)
      case array: WdlArray => WdlArray(array.wdlType, array.value map toUnixPath(docker))
      case map: WdlMap => WdlMap(map.wdlType, map.value mapValues toUnixPath(docker))
      case wdlValue => wdlValue
    }
  }

  def jobName: String = s"cromwell_${jobDescriptor.workflowDescriptor.id.shortString}_${jobDescriptor.call.unqualifiedName}"

  override def retryable = false

  lazy val workflowDescriptor = jobDescriptor.workflowDescriptor
  lazy val call = jobDescriptor.key.call
  lazy val pathBuilders = WorkflowPathsBackendInitializationData.pathBuilders(backendInitializationDataOption)
  lazy val callEngineFunction = SharedFileSystemExpressionFunctions(jobPaths, pathBuilders)
  override lazy val workflowId = jobDescriptor.workflowDescriptor.id
  override lazy val jobTag = jobDescriptor.key.tag

  lazy val isDockerRun = RuntimeAttributesValidation.extractOption(
    DockerValidation.instance, validatedRuntimeAttributes).isDefined

  def sharedReceive(jobOption: Option[SharedFileSystemJob]): Receive = LoggingReceive {
    case AbortJobCommand =>
      jobOption foreach tryKill
    case KvPutSuccess(_) => // expected after the KvPut in tellKvJobId
  }

  def instantiatedScript: String = {
    val pathTransformFunction = toUnixPath(isDockerRun) _
    val localizer = sharedFileSystem.localizeInputs(jobPaths.callInputsRoot, isDockerRun) _
    
    Command.instantiate(
      jobDescriptor, 
      callEngineFunction,
      localizer,
      pathTransformFunction
    ).get
  }

  override def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext) = {
    // Run now in receive, not in yet another Runnable.
    Future.fromTry(Try {
      mode match {
        case Execute =>
          tellMetadata(metadataKeyValues)
          executeScript()
        case Recover(recoveryId) =>
          recoveryId match {
            case job: SharedFileSystemJob => recoverScript(job)
            case other => throw new RuntimeException(s"Unable to recover $other")
          }
      }
    } recoverWith {
      case exception: Exception =>
        jobLogger.error(s"Error attempting to $mode the script", exception)
        Failure(exception)
    })
  }

  /**
    * Fire and forget info to the metadata service
    */
  def tellMetadata(metadataKeyValues: Map[String, Any]): Unit = {
    import cromwell.services.metadata.MetadataService.implicits.MetadataAutoPutter
    serviceRegistryActor.putMetadata(jobDescriptor.workflowDescriptor.id, Option(jobDescriptor.key), metadataKeyValues)
  }

  def executeScript(): ExecutionHandle = {
    val script = instantiatedScript
    jobLogger.info(s"`$script`")
    File(jobPaths.callExecutionRoot).createDirectories()
    val cwd = if (isDockerRun) jobPaths.callExecutionDockerRoot else jobPaths.callExecutionRoot
    writeScript(script, cwd)
    jobLogger.info(s"command: $processArgs")
    val runner = makeProcessRunner()
    val exitValue = runner.run()
    if (exitValue != 0) {
      FailedNonRetryableExecutionHandle(new RuntimeException("Unable to start job. " +
        s"Check the stderr file for possible errors: ${runner.stderrPath}"))
    } else {
      val runningJob = getJob(exitValue, runner.stdoutPath, runner.stderrPath)
      context.become(sharedReceive(Option(runningJob)) orElse super.receive)
      tellKvJobId(runningJob)
      jobLogger.info(s"job id: ${runningJob.jobId}")
      tellMetadata(Map("jobId" -> runningJob.jobId))
      SharedFileSystemPendingExecutionHandle(jobDescriptor, runningJob)
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
  private def writeScript(instantiatedCommand: String, cwd: Path) = {
    val rcPath = if (isDockerRun) jobPaths.toDockerPath(jobPaths.returnCode) else jobPaths.returnCode
    val rcTmpPath = s"$rcPath.tmp"

    val scriptBody = s"""

#!/bin/sh
(
 cd $cwd
 $instantiatedCommand
)
echo $$? > $rcTmpPath
mv $rcTmpPath $rcPath

""".trim + "\n"

    File(jobPaths.script).write(scriptBody)
  }

  /**
    * Send the job id of the running job to the key value store.
    *
    * @param runningJob The running job.
    */
  private def tellKvJobId(runningJob: SharedFileSystemJob): Unit = {
    val kvJobKey =
      KvJobKey(jobDescriptor.key.call.fullyQualifiedName, jobDescriptor.key.index, jobDescriptor.key.attempt)
    val scopedKey = ScopedKey(jobDescriptor.workflowDescriptor.id, kvJobKey, SharedFileSystemJob.JobIdKey)
    val kvValue = Option(runningJob.jobId)
    val kvPair = KvPair(scopedKey, kvValue)
    val kvPut = KvPut(kvPair)
    serviceRegistryActor ! kvPut
  }

  def recoverScript(job: SharedFileSystemJob): ExecutionHandle = {
    context.become(sharedReceive(Option(job)) orElse super.receive)
    // To avoid race conditions, check for the rc file after checking if the job is alive.
    if (isAlive(job) || File(jobPaths.returnCode).exists) {
      // If we're done, we'll get to the rc during the next poll.
      // Or if we're still running, return pending also.
      jobLogger.info(s"Recovering using job id: ${job.jobId}")
      SharedFileSystemPendingExecutionHandle(jobDescriptor, job)
    } else {
      // Could start executeScript(), but for now fail because we shouldn't be in this state.
      FailedNonRetryableExecutionHandle(new RuntimeException(
        s"Unable to determine that ${job.jobId} is alive, and ${jobPaths.returnCode} does not exist."), None)
    }
  }

  def isAlive(job: SharedFileSystemJob): Boolean = {
    val argv = checkAliveArgs(job).argv
    val stdout = pathPlusSuffix(jobPaths.stdout, "check")
    val stderr = pathPlusSuffix(jobPaths.stderr, "check")
    val checkAlive = new ProcessRunner(argv, stdout.path, stderr.path)
    checkAlive.run() == 0
  }

  def tryKill(job: SharedFileSystemJob): Unit = {
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

  def processReturnCode()(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    val returnCodeTry = Try(File(jobPaths.returnCode).contentAsString.stripLineEnd.toInt)

    lazy val badReturnCodeMessage = s"Call ${jobDescriptor.key.tag}: return code was ${returnCodeTry.getOrElse("(none)")}"

    lazy val badReturnCodeResponse = Future.successful(
      FailedNonRetryableExecutionHandle(new Exception(badReturnCodeMessage), returnCodeTry.toOption))

    lazy val abortResponse = Future.successful(AbortedExecutionHandle)

    def processSuccess(returnCode: Int) = {
      val successfulFuture = for {
        outputs <- Future.fromTry(processOutputs())
      } yield SuccessfulExecutionHandle(outputs, returnCode, jobPaths.detritusPaths, Seq.empty)

      successfulFuture recover {
        case failed: Throwable =>
          FailedNonRetryableExecutionHandle(failed, Option(returnCode))
      }
    }

    def stopFor(returnCode: Int) = {
      val continueOnReturnCode = RuntimeAttributesValidation.extract(
        ContinueOnReturnCodeValidation.instance, validatedRuntimeAttributes)
      !continueOnReturnCode.continueFor(returnCode)
    }

    def failForStderr = {
      val failOnStderr = RuntimeAttributesValidation.extract(
        FailOnStderrValidation.instance, validatedRuntimeAttributes)
      failOnStderr && File(jobPaths.stderr).size > 0
    }

    returnCodeTry match {
      case Success(SIGTERM) => abortResponse // Special case to check for SIGTERM exit code - implying abort
      case Success(SIGINT) => abortResponse // Special case to check for SIGINT exit code - implying abort
      case Success(returnCode) if stopFor(returnCode) => badReturnCodeResponse
      case Success(returnCode) if failForStderr => badReturnCodeResponse
      case Success(returnCode) => processSuccess(returnCode)
      case Failure(e) => badReturnCodeResponse
    }
  }

  override def poll(previous: ExecutionHandle)(implicit ec: ExecutionContext) = {
    previous match {
      case handle: SharedFileSystemPendingExecutionHandle =>
        val runId = handle.run
        jobLogger.debug(s"Polling Job $runId")
        File(jobPaths.returnCode).exists match {
          case true =>
            processReturnCode()
          case false =>
            jobLogger.debug(s"'${jobPaths.returnCode}' file does not exist yet")
            Future.successful(previous)
        }
      case failed: FailedNonRetryableExecutionHandle => Future.successful(failed)
      case successful: SuccessfulExecutionHandle => Future.successful(successful)
      case bad => Future.failed(new IllegalArgumentException(s"Unexpected execution handle: $bad"))
    }
  }

  private def processOutputs(): Try[JobOutputs] = {
    OutputEvaluator.evaluateOutputs(jobDescriptor, callEngineFunction, sharedFileSystem.outputMapper(jobPaths))
  }
}
