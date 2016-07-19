package cromwell.backend.sfs

import java.nio.file.{Path, Paths}

import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.event.LoggingReceive
import better.files._
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.async.AsyncBackendJobExecutionActor.ExecutionMode
import cromwell.backend.async.{AbortedExecutionHandle, AsyncBackendJobExecutionActor, ExecutionHandle, FailedNonRetryableExecutionHandle, NonRetryableExecution, SuccessfulExecutionHandle}
import cromwell.backend.io.{JobPaths, WorkflowPathsBackendInitializationData}
import cromwell.backend.validation._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptor, BackendJobDescriptorKey, ExecutionHash, OutputEvaluator}
import cromwell.core.PathFactory.EnhancedPath
import cromwell.core.logging.JobLogging
import cromwell.core.{JobOutputs, PathWriter}
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata._
import wdl4s.WdlExpression
import wdl4s.util.TryUtil
import wdl4s.values.{WdlArray, WdlFile, WdlMap, WdlValue}

import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.sys.process._
import scala.util.{Failure, Success, Try}

trait SharedFileSystemJob {
  def jobId: String
}

case class SharedFileSystemAsyncJobExecutionActorParams
(
  serviceRegistryActor: ActorRef,
  jobDescriptor: BackendJobDescriptor,
  configurationDescriptor: BackendConfigurationDescriptor,
  completionPromise: Promise[BackendJobExecutionResponse],
  supportsDocker: Boolean,
  runtimeAttributesBuilder: SharedFileSystemValidatedRuntimeAttributesBuilder,
  backendInitializationDataOption: Option[BackendInitializationData]
)

trait SharedFileSystemAsyncJobExecutionActor[JobType <: SharedFileSystemJob]
  extends Actor with ActorLogging with AsyncBackendJobExecutionActor with JobLogging {

  case class SharedFileSystemPendingExecutionHandle(jobDescriptor: BackendJobDescriptor,
                                                    run: JobType) extends ExecutionHandle {
    override val isDone = false
    override val result = NonRetryableExecution(new IllegalStateException(
      "SharedFileSystemPendingExecutionHandle cannot yield a result"))
  }

  context.become(sharedReceive(None) orElse super.receive)

  val SIGTERM = 143
  val SIGINT = 130

  override protected implicit def ec = context.dispatcher

  val params: SharedFileSystemAsyncJobExecutionActorParams

  def processArgs: SharedFileSystemCommand

  def getJob(exitValue: Int, stdout: Path, stderr: Path): JobType

  def killArgs(job: JobType): SharedFileSystemCommand

  def toPath(path: WdlValue): WdlValue = path

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

  def runsOnDocker = params.supportsDocker &&
    DockerValidation.optional.extract(validatedRuntimeAttributes).isDefined

  def jobName: String = s"cromwell_${jobDescriptor.descriptor.id.shortString}_${jobDescriptor.call.unqualifiedName}"

  def scriptDir: Path = if (runsOnDocker) jobPaths.callDockerRoot else jobPaths.callRoot

  override def retryable = false

  lazy val workflowDescriptor = jobDescriptor.descriptor
  lazy val call = jobDescriptor.key.call
  lazy val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobDescriptor.key)
  lazy val fileSystems = WorkflowPathsBackendInitializationData.fileSystems(backendInitializationDataOption)
  lazy val callEngineFunction = SharedFileSystemExpressionFunctions(jobPaths, fileSystems)
  override lazy val workflowId = jobDescriptor.descriptor.id
  override lazy val jobTag = jobDescriptor.key.tag

  private val lookup = jobDescriptor.inputs.apply _

  private def evaluate(wdlExpression: WdlExpression) = wdlExpression.evaluate(lookup, callEngineFunction)

  lazy val validatedRuntimeAttributes: ValidatedRuntimeAttributes = {
    val evaluateAttrs = call.task.runtimeAttributes.attrs mapValues evaluate
    // Fail the call if runtime attributes can't be evaluated
    val evaluatedAttributes = TryUtil.sequenceMap(evaluateAttrs, "Runtime attributes evaluation").get
    val builder = params.runtimeAttributesBuilder.withDockerSupport(params.supportsDocker)
    builder.build(evaluatedAttributes, jobDescriptor.descriptor.workflowOptions, jobLogger)
  }

  def sharedReceive(jobOption: Option[JobType]): Receive = LoggingReceive {
    case AbortJobCommand =>
      jobOption foreach tryKill
  }

  def instantiatedScript: String = {
    val pathTransformFunction: WdlValue => WdlValue = if (runsOnDocker) toDockerPath else toPath
    val tryCommand = sharedFileSystem.localizeInputs(jobPaths.callRoot,
      runsOnDocker, fileSystems, jobDescriptor.inputs) flatMap { localizedInputs =>
      call.task.instantiateCommand(localizedInputs, callEngineFunction, pathTransformFunction)
    }
    tryCommand.get
  }

  override def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext) = {
    Future.fromTry(Try {
      tellMetadata(startMetadataEvents)
      executeScript()
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

  def stdoutSubmit: PathWriter = pathPlusSuffix(jobPaths.stdout, "submit").untailed

  def stderrSubmit: PathWriter = pathPlusSuffix(jobPaths.stderr, "submit").untailed

  def executeScript(): SharedFileSystemPendingExecutionHandle = {
    val script = instantiatedScript
    jobLogger.info(s"`$script`")
    jobPaths.callRoot.createDirectories()
    writeScript(script, scriptDir.toAbsolutePath)
    jobLogger.info(s"command: $processArgs")
    val stdoutWriter = stdoutSubmit
    val stderrWriter = stderrSubmit
    val argv = processArgs.argv.map(_.toString)
    val proc = argv.run(ProcessLogger(stdoutWriter.writeWithNewline, stderrWriter.writeWithNewline))
    val run = getJob(proc, stdoutWriter, stderrWriter)
    context.become(sharedReceive(Option(run)) orElse super.receive)
    tellMetadata(Seq(metadataEvent("jobId", run.jobId)))
    SharedFileSystemPendingExecutionHandle(jobDescriptor, run)
  }

  def getJob(process: Process, stdoutWriter: PathWriter, stderrWriter: PathWriter): JobType = {
    val exitValue = process.exitValue()
    import cromwell.core.PathFactory.FlushingAndClosingWriter
    stdoutWriter.writer.flushAndClose()
    stderrWriter.writer.flushAndClose()
    getJob(exitValue, stdoutWriter.path, stderrWriter.path)
  }

  def tryKill(job: JobType): Unit = {
    import better.files._
    import cromwell.core.PathFactory._
    val returnCodeTmp = pathPlusSuffix(jobPaths.returnCode, "kill")
    val stdoutWriter = pathPlusSuffix(jobPaths.stdout, "kill").untailed
    val stderrWriter = pathPlusSuffix(jobPaths.stderr, "kill").untailed
    returnCodeTmp.write(s"$SIGTERM\n")
    returnCodeTmp.moveTo(jobPaths.returnCode)
    val argv = killArgs(job).argv.map(_.toString)
    val proc = argv.run(ProcessLogger(stdoutWriter.writeWithNewline, stderrWriter.writeWithNewline))
    proc.exitValue()
    stdoutWriter.writer.flushAndClose()
    stderrWriter.writer.flushAndClose()
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

    def stopFor(returnCode: Int) =
      !ContinueOnReturnCodeValidation.default.extract(validatedRuntimeAttributes).continueFor(returnCode)

    def failForStderr = FailOnStderrValidation.default.extract(validatedRuntimeAttributes) && jobPaths.stderr.size > 0

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
            cleanup()
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

  def cleanup(): Unit = ()

  /**
    * Writes the script file containing the user's command from the WDL as well
    * as some extra shell code for monitoring jobs
    */
  private def writeScript(instantiatedCommand: String, containerRoot: Path) = {
    jobPaths.script.write(
      s"""#!/bin/sh
          |cd $containerRoot
          |$instantiatedCommand
          |echo $$? > rc
          |""".stripMargin)
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
