package cromwell.backend.sfs

import java.nio.file.{FileAlreadyExistsException, Path, Paths}

import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.event.LoggingReceive
import better.files._
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.async.AsyncBackendJobExecutionActor.{ExecutionMode, JobId}
import cromwell.backend.async.{AbortedExecutionHandle, AsyncBackendJobExecutionActor, ExecutionHandle, FailedNonRetryableExecutionHandle, NonRetryableExecution, SuccessfulExecutionHandle}
import cromwell.backend.io.{JobPaths, WorkflowPathsBackendInitializationData}
import cromwell.backend.validation._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptor, BackendJobDescriptorKey, ExecutionHash, OutputEvaluator}
import cromwell.core.JobOutputs
import cromwell.core.logging.JobLogging
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata._
import wdl4s.WdlExpression
import wdl4s.util.TryUtil
import wdl4s.values.{WdlArray, WdlFile, WdlMap, WdlValue}

import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.util.{Failure, Success, Try}
import scala.concurrent.duration._

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
  extends Actor with ActorLogging with AsyncBackendJobExecutionActor with JobLogging {

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
    * Returns the command for killing a job.
    *
    * @param job The job to kill.
    * @return The command for killing a job.
    */
  def killArgs(job: SharedFileSystemJob): SharedFileSystemCommand

  override def jobDescriptor: BackendJobDescriptor = params.jobDescriptor

  override lazy val completionPromise: Promise[BackendJobExecutionResponse] = params.completionPromise

  lazy val serviceRegistryActor: ActorRef = params.serviceRegistryActor

  def configurationDescriptor: BackendConfigurationDescriptor = params.configurationDescriptor

  def backendInitializationDataOption: Option[BackendInitializationData] = params.backendInitializationDataOption

  def toDockerPath(path: WdlValue): WdlValue = {
    path match {
      case file: WdlFile => WdlFile(jobPaths.toDockerPath(Paths.get(path.valueString)).toAbsolutePath.toString)
      case array: WdlArray => WdlArray(array.wdlType, array.value map toDockerPath)
      case map: WdlMap => WdlMap(map.wdlType, map.value mapValues toDockerPath)
      case wdlValue => wdlValue
    }
  }

  def jobName: String = s"cromwell_${jobDescriptor.workflowDescriptor.id.shortString}_${jobDescriptor.call.unqualifiedName}"

  override def retryable = false

  lazy val workflowDescriptor = jobDescriptor.workflowDescriptor
  lazy val call = jobDescriptor.key.call
  lazy val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobDescriptor.key)
  lazy val fileSystems = WorkflowPathsBackendInitializationData.fileSystems(backendInitializationDataOption)
  lazy val callEngineFunction = SharedFileSystemExpressionFunctions(jobPaths, fileSystems)
  override lazy val workflowId = jobDescriptor.workflowDescriptor.id
  override lazy val jobTag = jobDescriptor.key.tag

  private val lookup = jobDescriptor.inputs.apply _

  private def evaluate(wdlExpression: WdlExpression) = wdlExpression.evaluate(lookup, callEngineFunction)

  private lazy val initializationData = BackendInitializationData.
    as[SharedFileSystemBackendInitializationData](backendInitializationDataOption)

  lazy val validatedRuntimeAttributes: ValidatedRuntimeAttributes = {
    val evaluateAttrs = call.task.runtimeAttributes.attrs mapValues evaluate
    // Fail the call if runtime attributes can't be evaluated
    val evaluatedAttributes = TryUtil.sequenceMap(evaluateAttrs, "Runtime attributes evaluation").get
    val builder = initializationData.runtimeAttributesBuilder
    builder.build(evaluatedAttributes, jobDescriptor.workflowDescriptor.workflowOptions, jobLogger)
  }

  lazy val isDockerRun = RuntimeAttributesValidation.extractOption(
    DockerValidation.instance, validatedRuntimeAttributes).isDefined

  def sharedReceive(jobOption: Option[SharedFileSystemJob]): Receive = LoggingReceive {
    case AbortJobCommand =>
      jobOption foreach tryKill
  }

  def instantiatedScript: String = {
    val pathTransformFunction: WdlValue => WdlValue = if (isDockerRun) toDockerPath else identity
    val tryCommand = sharedFileSystem.localizeInputs(jobPaths.callRoot,
      isDockerRun, fileSystems, jobDescriptor.inputs) flatMap { localizedInputs =>
      call.task.instantiateCommand(localizedInputs, callEngineFunction, pathTransformFunction)
    }
    tryCommand.get
  }

  override def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext) = {
    // Run now in receive, not in yet another Runnable.
    Future.fromTry(Try {
      tellMetadata(startMetadataEvents)
      executeScript()
    } recoverWith {
      case exception: Exception =>
        jobLogger.error("Error attempting to execute the script", exception)
        Failure(exception)
    })
  }

  private lazy val metadataJobKey = {
    val jobDescriptorKey: BackendJobDescriptorKey = jobDescriptor.key
    MetadataJobKey(jobDescriptorKey.call.fullyQualifiedName, jobDescriptorKey.index, jobDescriptorKey.attempt)
  }

  private def metadataKey(key: String) = MetadataKey(workflowId, Option(metadataJobKey), key)

  private def metadataEvent(key: String, value: Any) = MetadataEvent(metadataKey(key), MetadataValue(value))

  val runtimeAttributesEvents = validatedRuntimeAttributes.attributes map {
    case (key, value) =>
      metadataEvent(s"runtimeAttributes:$key", value)
  }

  def startMetadataEvents: Iterable[MetadataEvent] = runtimeAttributesEvents ++ List(
    metadataEvent(CallMetadataKeys.Stdout, jobPaths.stdout.toAbsolutePath),
    metadataEvent(CallMetadataKeys.Stderr, jobPaths.stderr.toAbsolutePath),
    // TODO: PBE: The REST endpoint toggles this value... how/where? Meanwhile, we read it decide to use the cache...
    metadataEvent("cache:allowResultReuse", true),
    metadataEvent(CallMetadataKeys.CallRoot, jobPaths.callRoot)
  )

  /**
    * Fire and forget info to the metadata service
    */
  def tellMetadata(events: Iterable[MetadataEvent]): Unit = {
    serviceRegistryActor ! PutMetadataAction(events)
  }

  def pathPlusSuffix(path: Path, suffix: String) = path.resolveSibling(s"${path.name}.$suffix")

  def executeScript(): ExecutionHandle = {
    val script = instantiatedScript
    jobLogger.info(s"`$script`")
    jobPaths.callRoot.createDirectories()
    val cwd = if (isDockerRun) jobPaths.callDockerRoot else jobPaths.callRoot
    writeScript(script, cwd)
    jobLogger.info(s"command: $processArgs")
    val runner = makeProcessRunner()
    val exitValue = runner.run()
    if (exitValue != 0) {
      FailedNonRetryableExecutionHandle(new RuntimeException("Unable to start job. " +
        s"Check the stderr file for possible errors: ${runner.stderrPath.fullPath}"))
    } else {
      val runningJob = getJob(exitValue, runner.stdoutPath, runner.stderrPath)
      context.become(sharedReceive(Option(runningJob)) orElse super.receive)
      jobLogger.info(s"job id: ${runningJob.jobId}")
      tellMetadata(Seq(metadataEvent("jobId", runningJob.jobId)))
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
    new ProcessRunner(processArgs.argv, jobPaths.stdout, jobPaths.stderr)
  }

  /**
    * Writes the script file containing the user's command from the WDL as well
    * as some extra shell code for monitoring jobs
    */
  private def writeScript(instantiatedCommand: String, cwd: Path) = {
    jobPaths.script.write(
      s"""#!/bin/sh
          |cd ${cwd.fullPath}
          |$instantiatedCommand
          |echo $$? > rc
          |""".stripMargin)
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
    val killer = new ProcessRunner(argv, stdout, stderr)
    killer.run()
  }

  def processReturnCode()(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    val returnCodeTry = Try(jobPaths.returnCode.contentAsString.stripLineEnd.toInt)

    lazy val badReturnCodeMessage =
      s"Call ${call.fullyQualifiedName}, Workflow ${workflowDescriptor.id}: " +
        s"return code was ${returnCodeTry.getOrElse("(none)")}"

    lazy val badReturnCodeResponse = Future.successful(
      FailedNonRetryableExecutionHandle(new Exception(badReturnCodeMessage), returnCodeTry.toOption))

    lazy val abortResponse = Future.successful(AbortedExecutionHandle)

    def processSuccess(returnCode: Int) = {
      val successfulFuture = for {
        outputs <- Future.fromTry(processOutputs())
        hash <- ExecutionHash.completelyRandomExecutionHash
      } yield SuccessfulExecutionHandle(outputs, returnCode, hash, None)

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
      failOnStderr && jobPaths.stderr.size > 0
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
        jobPaths.returnCode.exists match {
          case true =>
            processReturnCode()
          case false =>
            jobLogger.info(s"'${jobPaths.returnCode}' file does not exist yet")
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

  private val sharedFileSystem = new SharedFileSystem {
    override lazy val sharedFileSystemConfig = {
      import lenthall.config.ScalaConfig._
      val config = configurationDescriptor.backendConfig.getConfigOption("filesystems.local")
      config.getOrElse(ConfigFactory.parseString("""localization: [ "hard-link", "soft-link", "copy" ]"""))
    }
  }
}
