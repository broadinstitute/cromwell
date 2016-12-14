package cromwell.backend.impl.htcondor

import java.nio.file.attribute.PosixFilePermission
import java.util.UUID

import akka.actor.{ActorRef, Props}
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, JobFailedNonRetryableResponse, JobSucceededResponse}
import cromwell.backend._
import cromwell.backend.impl.htcondor.caching.CacheActor._
import cromwell.backend.impl.htcondor.caching.localization.CachedResultLocalization
import cromwell.backend.io.JobPathsWithDocker
import cromwell.backend.sfs.{SharedFileSystem, SharedFileSystemExpressionFunctions}
import cromwell.backend.wdl.Command
import cromwell.core.path.JavaWriterImplicits._
import cromwell.core.path.{DefaultPathBuilder, PathBuilder}
import cromwell.services.keyvalue.KeyValueServiceActor._
import cromwell.services.metadata.CallMetadataKeys
import lenthall.util.TryUtil
import org.apache.commons.codec.digest.DigestUtils
import wdl4s.EvaluatedTaskInputs
import wdl4s.parser.MemoryUnit
import wdl4s.types.{WdlArrayType, WdlFileType}
import wdl4s.values.WdlArray

import scala.concurrent.{Future, Promise}
import scala.sys.process.ProcessLogger
import scala.util.{Failure, Success, Try}

object HtCondorJobExecutionActor {
  val HtCondorJobIdKey = "htCondor_job_id"

  val pathBuilders = List(DefaultPathBuilder)

  def props(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor, serviceRegistryActor: ActorRef, cacheActorProps: Option[Props]): Props =
    Props(new HtCondorJobExecutionActor(jobDescriptor, configurationDescriptor, serviceRegistryActor, cacheActorProps))
}

class HtCondorJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                                override val configurationDescriptor: BackendConfigurationDescriptor,
                                serviceRegistryActor: ActorRef,
                                cacheActorProps: Option[Props])
  extends BackendJobExecutionActor with CachedResultLocalization with SharedFileSystem {

  import HtCondorJobExecutionActor._
  import better.files._

  private val tag = s"CondorJobExecutionActor-${jobDescriptor.call.fullyQualifiedName}:"
  override val pathBuilders: List[PathBuilder] = HtCondorJobExecutionActor.pathBuilders

  implicit val executionContext = context.dispatcher

  lazy val cmds = new HtCondorCommands
  lazy val extProcess = new HtCondorProcess

  private val fileSystemsConfig = configurationDescriptor.backendConfig.getConfig("filesystems")
  override val sharedFileSystemConfig = fileSystemsConfig.getConfig("local")
  private val workflowDescriptor = jobDescriptor.workflowDescriptor
  private val jobPaths = new JobPathsWithDocker(jobDescriptor.key, workflowDescriptor, configurationDescriptor.backendConfig)

  // Files
  private val executionDir = jobPaths.callExecutionRoot
  private val returnCodePath = jobPaths.returnCode
  private val stdoutPath = jobPaths.stdout
  private val stderrPath = jobPaths.stderr
  private val scriptPath = jobPaths.script

  // stdout stderr writers for submit file logs
  private val submitFilePath = executionDir.resolve("submitfile")
  private val submitFileStderr = executionDir.resolve("submitfile.stderr")
  private val submitFileStdout = executionDir.resolve("submitfile.stdout")
  private val htCondorLog = executionDir.resolve(s"${jobDescriptor.call.unqualifiedName}.log")

  private lazy val stdoutWriter = extProcess.untailedWriter(submitFileStdout)
  private lazy val stderrWriter = extProcess.tailedWriter(100, submitFileStderr)

  private val call = jobDescriptor.key.call
  private val callEngineFunction = SharedFileSystemExpressionFunctions(jobPaths, pathBuilders)

  private val lookup = jobDescriptor.fullyQualifiedInputs.apply _

  private val runtimeAttributes = {
    val evaluateAttrs = call.task.runtimeAttributes.attrs mapValues (_.evaluate(lookup, callEngineFunction))
    // Fail the call if runtime attributes can't be evaluated
    val runtimeMap = TryUtil.sequenceMap(evaluateAttrs, "Runtime attributes evaluation").get
    HtCondorRuntimeAttributes(runtimeMap, jobDescriptor.workflowDescriptor.workflowOptions)
  }

  private val cacheActor = cacheActorProps match {
    case Some(props) => Some(context.actorOf(props, s"CacheActor-${jobDescriptor.call.fullyQualifiedName}"))
    case None => None
  }

  log.debug("{} Calculating hash for current job.", tag)
  lazy private val jobHash = calculateHash

  private val executionResponse = Promise[BackendJobExecutionResponse]()

  // Message sent (by self, to self) wrapping over the response produced by HtCondor
  private final case class JobExecutionResponse(resp: BackendJobExecutionResponse)

  // Message sent (by self, to self) to trigger a status check to HtCondor
  private final case class TrackTaskStatus(id: String)

  private var condorJobId: Option[String] = None

  private val pollingInterval = configurationDescriptor.backendConfig.getInt("poll-interval")

  override def receive = super.receive orElse {
    case JobExecutionResponse(resp) =>
      log.debug("{}: Completing job [{}] with response: [{}]", tag, jobDescriptor.key, resp)
      executionResponse trySuccess resp
      ()
    case TrackTaskStatus(id) =>
      // Avoid the redundant status check if the response is already completed (e.g. in case of abort)
      if (!executionResponse.isCompleted) trackTask(id)

    // Messages received from Caching actor
    case ExecutionResultFound(succeededResponse) =>
      executionResponse trySuccess localizeCachedResponse(succeededResponse)
      ()
    case ExecutionResultNotFound => prepareAndExecute()
    case ExecutionResultStored(hash) => log.debug("{} Cache entry was stored for Job with hash {}.", tag, hash)
    case ExecutionResultAlreadyExist => log.warning("{} Cache entry for hash {} already exist.", tag, jobHash)

    // Messages received from KV actor
    case KvPair(scopedKey, Some(jobId)) if scopedKey.key == HtCondorJobIdKey =>
      log.info("{} Found job id {}. Trying to recover job now.", tag, jobId)
      self ! TrackTaskStatus(jobId)
    case KvKeyLookupFailed(_) =>
      log.debug("{} Job id not found. Falling back to execute.", tag)
      execute
      // -Ywarn-value-discard
      ()
    case KvFailure(_, e) =>
      log.error("{} Failure attempting to look up HtCondor job id. Exception message: {}. Falling back to execute.", tag, e.getMessage)
      execute
      // -Ywarn-value-discard
      ()
  }

  /**
    * Restart or resume a previously-started job.
    */
  override def recover: Future[BackendJobExecutionResponse] = {
    log.warning("{} Trying to recover job {}.", tag, jobDescriptor.key.call.fullyQualifiedName)
    serviceRegistryActor ! KvGet(ScopedKey(jobDescriptor.workflowDescriptor.id,
      KvJobKey(jobDescriptor.key.call.fullyQualifiedName, jobDescriptor.key.index, jobDescriptor.key.attempt),
      HtCondorJobIdKey))
    executionResponse.future
  }

  /**
    * Execute a new job.
    */
  override def execute: Future[BackendJobExecutionResponse] = {
    log.debug("{} Checking if hash {{}} is in the cache.", tag, jobHash)
    cacheActor match {
      case Some(actorRef) => actorRef ! ReadExecutionResult(jobHash)
      case None => prepareAndExecute()
    }
    executionResponse.future
  }

  /**
    * Abort a running job.
    */
  override def abort(): Unit = {
    // Nothing to do in case `condorJobId` is not defined
    condorJobId foreach { id =>
      log.info("{}: Aborting job [{}:{}].", tag, jobDescriptor.key.tag, id)
      val abortProcess = new HtCondorProcess
      val argv = Seq(HtCondorCommands.Remove, id)
      val process = abortProcess.externalProcess(argv)
      val exitVal = process.exitValue()
      if (exitVal == 0)
        log.info("{}: Job {} successfully killed and removed from the queue.", tag, id)
      else
        log.error("{}: Failed to kill / remove job {}. Exit Code: {}, Stderr: {}", tag, id, exitVal, abortProcess.processStderr)
    }
  }

  private def executeTask(): Unit = {
    val argv = Seq(HtCondorCommands.Submit, submitFilePath.toString)
    val process = extProcess.externalProcess(argv, ProcessLogger(stdoutWriter.writeWithNewline, stderrWriter.writeWithNewline))
    val condorReturnCode = process.exitValue() // blocks until process (i.e. condor submission) finishes
    log.debug("{} Return code of condor submit command: {}", tag, condorReturnCode)

    List(stdoutWriter.writer, stderrWriter.writer).foreach(_.flushAndClose())

    condorReturnCode match {
      case 0 if File(submitFileStderr).lines.toList.isEmpty =>
        log.info("{} {} submitted to HtCondor. Waiting for the job to complete via. RC file status.", tag, jobDescriptor.call.fullyQualifiedName)
        val job = HtCondorCommands.SubmitOutputPattern.r
        //Number of lines in stdout for submit job will be 3 at max therefore reading all lines at once.
        log.debug(s"{} Output of submit process : {}", tag, File(submitFileStdout).lines.toList)
        val line = File(submitFileStdout).lines.toList.last
        line match {
          case job(jobId, clusterId) =>
            val overallJobIdentifier = s"$clusterId.${jobId.toInt - 1}" // Condor has 0 based indexing on the jobs, probably won't work on stuff like `queue 150`
            log.info("{} {} mapped to HtCondor JobID: {}", tag, jobDescriptor.call.fullyQualifiedName, overallJobIdentifier)
            serviceRegistryActor ! KvPut(KvPair(ScopedKey(jobDescriptor.workflowDescriptor.id,
              KvJobKey(jobDescriptor.key.call.fullyQualifiedName, jobDescriptor.key.index, jobDescriptor.key.attempt),
              HtCondorJobIdKey), Option(overallJobIdentifier)))
            condorJobId = Option(overallJobIdentifier)
            self ! TrackTaskStatus(overallJobIdentifier)

          case _ => self ! JobExecutionResponse(JobFailedNonRetryableResponse(jobDescriptor.key,
            new IllegalStateException("Failed to retrieve job(id) and cluster id"), Option(condorReturnCode)))
        }

      case 0 =>
        log.error(s"Unexpected! Received return code for condor submission as 0, although stderr file is non-empty: {}", File(submitFileStderr).lines)
        self ! JobExecutionResponse(JobFailedNonRetryableResponse(jobDescriptor.key,
          new IllegalStateException(s"Execution process failed. HtCondor returned zero status code but non empty stderr file: $condorReturnCode"),
          Option(condorReturnCode)))

      case nonZeroExitCode: Int =>
        self ! JobExecutionResponse(JobFailedNonRetryableResponse(jobDescriptor.key,
          new IllegalStateException(s"Execution process failed. HtCondor returned non zero status code: $condorReturnCode"), Option(condorReturnCode)))
    }
  }

  private def trackTask(jobIdentifier: String): Unit = {
    val jobReturnCode = Try(extProcess.jobReturnCode(jobIdentifier, returnCodePath))
    log.debug("{} Process complete. RC file now exists with value: {}", tag, jobReturnCode)

    jobReturnCode match {
      case Success(None) =>
        import scala.concurrent.duration._
        // Job is still running in HtCondor. Check back again after `pollingInterval` seconds
        context.system.scheduler.scheduleOnce(pollingInterval.seconds, self, TrackTaskStatus(jobIdentifier))
        ()
      case Success(Some(rc)) if runtimeAttributes.continueOnReturnCode.continueFor(rc) => self ! JobExecutionResponse(processSuccess(rc))
      case Success(Some(rc)) => self ! JobExecutionResponse(JobFailedNonRetryableResponse(jobDescriptor.key,
        new IllegalStateException("Job exited with invalid return code: " + rc), Option(rc)))
      case Failure(error) => self ! JobExecutionResponse(JobFailedNonRetryableResponse(jobDescriptor.key, error, None))
    }
  }

  private def processSuccess(rc: Int): BackendJobExecutionResponse = {
    evaluateOutputs(callEngineFunction, outputMapper(jobPaths)) match {
      case Success(outputs) =>
        val succeededResponse = JobSucceededResponse(jobDescriptor.key, Some(rc), outputs, None, Seq.empty)
        log.debug("{} Storing data into cache for hash {}.", tag, jobHash)
        // If cache fails to store data for any reason it should not stop the workflow/task execution but log the issue.
        cacheActor foreach { _ ! StoreExecutionResult(jobHash, succeededResponse) }
        succeededResponse
      case Failure(e) =>
        val message = Option(e.getMessage) map {
          ": " + _
        } getOrElse ""
        JobFailedNonRetryableResponse(jobDescriptor.key, new Throwable("Failed post processing of outputs" + message, e), Option(rc))
    }
  }

  private def calculateHash: String = {
    val cmd = Command.instantiate(jobDescriptor, callEngineFunction) match {
      case Success(command) => command
      case Failure(ex) =>
        val errMsg = s"$tag Cannot instantiate job command for caching purposes due to ${ex.getMessage}."
        log.error(ex.getCause, errMsg)
        throw new IllegalStateException(errMsg, ex.getCause)
    }
    val str = Seq[Any](cmd,
      runtimeAttributes.failOnStderr,
      runtimeAttributes.dockerImage.getOrElse(""),
      runtimeAttributes.dockerWorkingDir.getOrElse(""),
      runtimeAttributes.dockerOutputDir.getOrElse(""),
      runtimeAttributes.cpu.toString,
      runtimeAttributes.memory.toString,
      runtimeAttributes.disk.toString).mkString
    DigestUtils.md5Hex(str)
  }

  private def createExecutionFolderAndScript(): Unit = {
    try {
      log.debug("{} Creating execution folder: {}", tag, executionDir)
      executionDir.toString.toFile.createIfNotExists(asDirectory = true, createParents = true)

      log.debug("{} Resolving job command", tag)
      val command = localizeInputs(jobPaths.callInputsRoot, runtimeAttributes.dockerImage.isDefined)(jobDescriptor.inputDeclarations) flatMap {
        localizedInputs => resolveJobCommand(localizedInputs)
      }

      log.debug("{} Creating bash script for executing command: {}", tag, command)
      cmds.writeScript(command.get, scriptPath.toAbsolutePath, executionDir.toAbsolutePath) // Writes the bash script for executing the command
      File(scriptPath).addPermission(PosixFilePermission.OWNER_EXECUTE) // Add executable permissions to the script.
      //TODO: Need to append other runtime attributes from Wdl to Condor submit file
      val attributes: Map[String, Any] = Map(HtCondorRuntimeKeys.Executable -> scriptPath.toAbsolutePath,
          HtCondorRuntimeKeys.InitialWorkingDir -> jobPaths.callExecutionRoot.toAbsolutePath,
          HtCondorRuntimeKeys.Output -> stdoutPath.toAbsolutePath,
          HtCondorRuntimeKeys.Error -> stderrPath.toAbsolutePath,
          HtCondorRuntimeKeys.Log -> htCondorLog.toAbsolutePath,
          HtCondorRuntimeKeys.LogXml -> true,
          HtCondorRuntimeKeys.LeaveInQueue -> true,
          HtCondorRuntimeKeys.Cpu -> runtimeAttributes.cpu,
          HtCondorRuntimeKeys.Memory -> runtimeAttributes.memory.to(MemoryUnit.MB).amount.toLong,
          HtCondorRuntimeKeys.Disk -> runtimeAttributes.disk.to(MemoryUnit.KB).amount.toLong
        )

      cmds.generateSubmitFile(submitFilePath, attributes, runtimeAttributes.nativeSpecs) // This writes the condor submit file
      ()

    } catch {
      case ex: Exception =>
        log.error(ex, "Failed to prepare task: " + ex.getMessage)
        throw ex
    }
  }

  private def resolveJobCommand(localizedInputs: EvaluatedTaskInputs): Try[String] = {
    val command = if (runtimeAttributes.dockerImage.isDefined) {
      modifyCommandForDocker(call.task.instantiateCommand(localizedInputs, callEngineFunction, identity), localizedInputs)
    } else {
      call.task.instantiateCommand(localizedInputs, callEngineFunction, identity)
    }
    command match {
      case Success(cmd) => tellMetadata(Map("command" -> cmd))
      case Failure(ex) =>
        log.error("{} failed to resolve command due to exception:{}", tag, ex)
        tellMetadata(Map(s"${CallMetadataKeys.Failures}[${UUID.randomUUID().toString}]" -> ex.getMessage))
    }
    command
  }

  /**
    * Fire and forget data to the metadata service
    */
  private def tellMetadata(metadataKeyValues: Map[String, Any]): Unit = {
    import cromwell.services.metadata.MetadataService.implicits.MetadataAutoPutter
    serviceRegistryActor.putMetadata(jobDescriptor.workflowDescriptor.id, Option(jobDescriptor.key), metadataKeyValues)
  }

  private def modifyCommandForDocker(jobCmd: Try[String], localizedInputs: EvaluatedTaskInputs): Try[String] = {
    Try {
      val dockerInputDataVol = localizedInputs.values.collect {
        case file if file.wdlType == WdlFileType =>
          val limit = file.valueString.lastIndexOf("/")
          Seq(file.valueString.substring(0, limit))
        case files if files.wdlType == WdlArrayType(WdlFileType) => files.asInstanceOf[WdlArray].value map { file =>
          val limit = file.valueString.lastIndexOf("/")
          file.valueString.substring(0, limit)
        }
      }.flatten.toSeq

      log.debug("{} List of input volumes: {}", tag, dockerInputDataVol.mkString(","))
      val dockerCmd = configurationDescriptor.backendConfig.getString("docker.cmd")
      val defaultWorkingDir = configurationDescriptor.backendConfig.getString("docker.defaultWorkingDir")
      val defaultOutputDir = configurationDescriptor.backendConfig.getString("docker.defaultOutputDir")
      val dockerVolume = "-v %s:%s"
      val dockerVolumeInputs = s"$dockerVolume:ro"
      // `v.get` is safe below since we filtered the list earlier with only defined elements
      val inputVolumes = dockerInputDataVol.distinct.map(v => dockerVolumeInputs.format(v, v)).mkString(" ")
      val outputVolume = dockerVolume.format(executionDir.toAbsolutePath.toString, runtimeAttributes.dockerOutputDir.getOrElse(defaultOutputDir))
      val workingDir = dockerVolume.format(executionDir.toAbsolutePath.toString, runtimeAttributes.dockerWorkingDir.getOrElse(defaultWorkingDir))
      val cmd = dockerCmd.format(runtimeAttributes.dockerWorkingDir.getOrElse(defaultWorkingDir), workingDir, inputVolumes, outputVolume, runtimeAttributes.dockerImage.get, jobCmd.get)
      log.debug("{} Docker command line to be used for task execution: {}.", tag, cmd)
      cmd
    }
  }

  private def prepareAndExecute(): Unit = {
    try {
      createExecutionFolderAndScript()
      executeTask()
    } catch {
      case e: Exception => self ! JobExecutionResponse(JobFailedNonRetryableResponse(jobDescriptor.key, e, None))
    }
  }

  private def localizeCachedResponse(succeededResponse: JobSucceededResponse): BackendJobExecutionResponse = {
    Try(localizeCachedOutputs(executionDir, succeededResponse.jobOutputs)) match {
      case Success(outputs) =>
        executionDir.toString.toFile.createIfNotExists(asDirectory = true, createParents = true)
        JobSucceededResponse(jobDescriptor.key, succeededResponse.returnCode, outputs, None, Seq.empty)
      case Failure(exception) => JobFailedNonRetryableResponse(jobDescriptor.key, exception, None)
    }
  }
}
