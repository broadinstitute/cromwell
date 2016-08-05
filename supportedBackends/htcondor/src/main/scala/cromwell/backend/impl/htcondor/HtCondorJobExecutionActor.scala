package cromwell.backend.impl.htcondor

import java.nio.file.attribute.PosixFilePermission
import java.nio.file.FileSystems

import akka.actor.{ActorRef, Props}
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, FailedNonRetryableResponse, SucceededResponse}
import cromwell.backend._
import cromwell.backend.impl.htcondor.caching.CacheActor._
import cromwell.backend.io.JobPaths
import cromwell.backend.sfs.{SharedFileSystem, SharedFileSystemExpressionFunctions}
import cromwell.services.keyvalue.KeyValueServiceActor._
import org.apache.commons.codec.digest.DigestUtils
import wdl4s._
import wdl4s.parser.MemoryUnit
import wdl4s.types.{WdlArrayType, WdlFileType}
import wdl4s.util.TryUtil
import wdl4s.values.WdlArray

import scala.concurrent.{Future, Promise}
import scala.sys.process.ProcessLogger
import scala.util.{Failure, Success, Try}
import scala.language.postfixOps

object HtCondorJobExecutionActor {
  val HtCondorJobIdKey = "htCondor_job_id"

  val fileSystems = List(FileSystems.getDefault)

  def props(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor, serviceRegistryActor: ActorRef, cacheActorProps: Option[Props]): Props =
    Props(new HtCondorJobExecutionActor(jobDescriptor, configurationDescriptor, serviceRegistryActor, cacheActorProps))
}

class HtCondorJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                                override val configurationDescriptor: BackendConfigurationDescriptor,
                                serviceRegistryActor: ActorRef,
                                cacheActorProps: Option[Props]) extends BackendJobExecutionActor with SharedFileSystem {

  import HtCondorJobExecutionActor._
  import better.files._
  import cromwell.core.PathFactory._

  private val tag = s"CondorJobExecutionActor-${jobDescriptor.call.fullyQualifiedName}:"

  lazy val cmds = new HtCondorCommands
  lazy val extProcess = new HtCondorProcess

  private val fileSystemsConfig = configurationDescriptor.backendConfig.getConfig("filesystems")
  override val sharedFileSystemConfig = fileSystemsConfig.getConfig("local")
  private val workflowDescriptor = jobDescriptor.descriptor
  private val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobDescriptor.key)

  // Files
  private val executionDir = jobPaths.callRoot
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
  private val callEngineFunction = SharedFileSystemExpressionFunctions(jobPaths, fileSystems)

  private val lookup = jobDescriptor.inputs.apply _

  private val runtimeAttributes = {
    val evaluateAttrs = call.task.runtimeAttributes.attrs mapValues (_.evaluate(lookup, callEngineFunction))
    // Fail the call if runtime attributes can't be evaluated
    val runtimeMap = TryUtil.sequenceMap(evaluateAttrs, "Runtime attributes evaluation").get
    HtCondorRuntimeAttributes(runtimeMap, jobDescriptor.descriptor.workflowOptions)
  }

  private val cacheActor = cacheActorProps match {
    case Some(props) => Some(context.actorOf(props, s"CacheActor-${jobDescriptor.call.fullyQualifiedName}"))
    case None => None
  }

  log.debug("{} Calculating hash for current job.", tag)
  lazy private val jobHash = calculateHash

  private val executionResponse = Promise[BackendJobExecutionResponse]()

  override def receive = super.receive orElse {
    case ExecutionResultFound(succeededResponse) => executionResponse success succeededResponse.copy(jobKey = jobDescriptor.key)

    case ExecutionResultNotFound => prepareAndExecute

    case ExecutionResultStored(hash) => log.debug("{} Cache entry was stored for Job with hash {}.", tag, hash)

    case ExecutionResultAlreadyExist => log.warning("{} Cache entry for hash {} already exist.", tag, jobHash)

    case KvPair(scopedKey, Some(jobId)) if scopedKey.key == HtCondorJobIdKey =>
      log.info("{} Found job id {}. Trying to recover job now.", tag, jobId)
      executionResponse success trackTaskToCompletion(jobId)

    case KvKeyLookupFailed(_) =>
      log.debug("{} Job id not found. Falling back to execute.", tag)
      execute

    case KvFailure(_, e) =>
      log.error(s"$tag Failure attempting to look up HtCondor job id. Exception message: ${e.getMessage}. Falling back to execute.")
      execute
  }

  /**
    * Restart or resume a previously-started job.
    */
  override def recover: Future[BackendJobExecutionResponse] = {
    log.warning("{} Trying to recover job {}.", tag, jobDescriptor.key.call.fullyQualifiedName)
    serviceRegistryActor ! KvGet(ScopedKey(jobDescriptor.descriptor.id, jobDescriptor.key, HtCondorJobIdKey))
    executionResponse.future
  }

  /**
    * Execute a new job.
    */
  override def execute: Future[BackendJobExecutionResponse] = {
    log.debug("{} Checking if hash {{}} is in the cache.", tag, jobHash)
    cacheActor match {
      case Some(actorRef) => actorRef ! ReadExecutionResult(jobHash)
      case None => prepareAndExecute
    }
    executionResponse.future
  }

  /**
    * Abort a running job.
    */
  override def abort(): Unit = throw new UnsupportedOperationException("HtCondorBackend currently doesn't support aborting jobs.")

  private def executeTask(): BackendJobExecutionResponse = {
    val argv = Seq(HtCondorCommands.Submit, submitFilePath.toString)
    val process = extProcess.externalProcess(argv, ProcessLogger(stdoutWriter.writeWithNewline, stderrWriter.writeWithNewline))
    val condorReturnCode = process.exitValue() // blocks until process (i.e. condor submission) finishes
    log.debug("{} Return code of condor submit command: {}", tag, condorReturnCode)

    List(stdoutWriter.writer, stderrWriter.writer).foreach(_.flushAndClose())

    condorReturnCode match {
      case 0 if submitFileStderr.lines.toList.isEmpty =>
        log.info("{} {} submitted to HtCondor. Waiting for the job to complete via. RC file status.", tag, jobDescriptor.call.fullyQualifiedName)
        val job = HtCondorCommands.SubmitOutputPattern.r
        //Number of lines in stdout for submit job will be 3 at max therefore reading all lines at once.
        log.debug(s"{} Output of submit process : {}", tag, submitFileStdout.lines.toList)
        val line = submitFileStdout.lines.toList.last
        line match {
          case job(jobId, clusterId) =>
            val overallJobIdentifier = s"$clusterId.${jobId.toInt - 1}" // Condor has 0 based indexing on the jobs, probably won't work on stuff like `queue 150`
            log.info("{} {} mapped to HtCondor JobID: {}", tag, jobDescriptor.call.fullyQualifiedName, overallJobIdentifier)
            serviceRegistryActor ! KvPut(KvPair(ScopedKey(jobDescriptor.descriptor.id, jobDescriptor.key, HtCondorJobIdKey), Option(overallJobIdentifier)))
            trackTaskToCompletion(overallJobIdentifier)

          case _ => FailedNonRetryableResponse(jobDescriptor.key,
            new IllegalStateException("Failed to retrieve job(id) and cluster id"), Option(condorReturnCode))
        }

      case 0 =>
        log.error(s"Unexpected! Received return code for condor submission as 0, although stderr file is non-empty: {}", submitFileStderr.lines)
        FailedNonRetryableResponse(jobDescriptor.key,
          new IllegalStateException(s"Execution process failed. HtCondor returned zero status code but non empty stderr file: $condorReturnCode"), Option(condorReturnCode))

      case nonZeroExitCode: Int =>
        FailedNonRetryableResponse(jobDescriptor.key,
          new IllegalStateException(s"Execution process failed. HtCondor returned non zero status code: $condorReturnCode"), Option(condorReturnCode))
    }
  }

  private def trackTaskToCompletion(jobIdentifier: String): BackendJobExecutionResponse = {
    val jobReturnCode = Try(extProcess.jobReturnCode(jobIdentifier, returnCodePath)) // Blocks until process completes
    log.debug("{} Process complete. RC file now exists with value: {}", tag, jobReturnCode)

    // TODO: Besides return code, do we also need to check based on stderr file?
    jobReturnCode match {
      case Success(rc) if rc == 0 | runtimeAttributes.continueOnReturnCode.continueFor(rc) => processSuccess(rc)
      case Success(rc) => FailedNonRetryableResponse(jobDescriptor.key,
        new IllegalStateException("Job exited with invalid return code: " + rc), Option(rc))
      case Failure(error) => FailedNonRetryableResponse(jobDescriptor.key, error, None)
    }
  }

  private def processSuccess(rc: Int) = {
    evaluateOutputs(callEngineFunction, outputMapper(jobPaths)) match {
      case Success(outputs) =>
        val succeededResponse = SucceededResponse(jobDescriptor.key, Some(rc), outputs)
        log.debug("{} Storing data into cache for hash {}.", tag, jobHash)
        //If cache fails to store data for any reason it should not stop the workflow/task execution but log the issue.
        cacheActor foreach { _ ! StoreExecutionResult(jobHash, succeededResponse) }
        succeededResponse
      case Failure(e) =>
        val message = Option(e.getMessage) map {
          ": " + _
        } getOrElse ""
        FailedNonRetryableResponse(jobDescriptor.key, new Throwable("Failed post processing of outputs" + message, e), Option(rc))
    }
  }

  private def calculateHash: String = {
    val cmd = call.task.instantiateCommand(jobDescriptor.inputs, callEngineFunction, identity) match {
      case Success(cmd) => cmd
      case Failure(ex) =>
        val errMsg = s"$tag Cannot instantiate job command for caching purposes due to ${ex.getMessage}."
        log.error(ex.getCause, errMsg)
        throw new IllegalStateException(errMsg, ex.getCause)
    }
    val str = Seq(cmd,
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
      executionDir.toString.toFile.createIfNotExists(true)

      log.debug("{} Resolving job command", tag)
      val command = localizeInputs(jobPaths.callRoot, runtimeAttributes.dockerImage.isDefined, fileSystems, jobDescriptor.inputs) flatMap {
        localizedInputs => resolveJobCommand(localizedInputs)
      }

      log.debug("{} Creating bash script for executing command: {}", tag, command)
      cmds.writeScript(command.get, scriptPath.toAbsolutePath, executionDir.toAbsolutePath) // Writes the bash script for executing the command
      scriptPath.addPermission(PosixFilePermission.OWNER_EXECUTE) // Add executable permissions to the script.
      //TODO: Need to append other runtime attributes from Wdl to Condor submit file
      val attributes: Map[String, Any] = Map(HtCondorRuntimeKeys.Executable -> scriptPath.toAbsolutePath,
          HtCondorRuntimeKeys.InitialWorkingDir -> jobPaths.callRoot.toAbsolutePath,
          HtCondorRuntimeKeys.Output -> stdoutPath.toAbsolutePath,
          HtCondorRuntimeKeys.Error -> stderrPath.toAbsolutePath,
          HtCondorRuntimeKeys.Log -> htCondorLog.toAbsolutePath,
          HtCondorRuntimeKeys.LogXml -> true,
          HtCondorRuntimeKeys.LeaveInQueue -> true,
          HtCondorRuntimeKeys.Cpu -> runtimeAttributes.cpu,
          HtCondorRuntimeKeys.Memory -> runtimeAttributes.memory.to(MemoryUnit.MB).amount.toLong,
          HtCondorRuntimeKeys.Disk -> runtimeAttributes.disk.to(MemoryUnit.KB).amount.toLong
        )

      cmds.generateSubmitFile(submitFilePath, attributes) // This writes the condor submit file

    } catch {
      case ex: Exception =>
        log.error(ex, "Failed to prepare task: " + ex.getMessage)
        throw ex
    }
  }

  private def resolveJobCommand(localizedInputs: CallInputs): Try[String] = {
    if (runtimeAttributes.dockerImage.isDefined)
      modifyCommandForDocker(call.task.instantiateCommand(localizedInputs, callEngineFunction, identity), localizedInputs)
    else
      call.task.instantiateCommand(localizedInputs, callEngineFunction, identity)
  }

  private def modifyCommandForDocker(jobCmd: Try[String], localizedInputs: CallInputs): Try[String] = {
    Try {
      val dockerInputDataVol = localizedInputs.collect {
        case (k, v) => v match {
          case file if file.wdlType == WdlFileType =>
            val limit = file.valueString.lastIndexOf("/")
            Seq(file.valueString.substring(0, limit))
          case files if files.wdlType == WdlArrayType(WdlFileType) => files.asInstanceOf[WdlArray].value map { file =>
            val limit = file.valueString.lastIndexOf("/")
            file.valueString.substring(0, limit)
          }
        }
      }.flatten.toSeq

      log.debug("{} List of input volumes: {}", tag, dockerInputDataVol.mkString(","))
      val dockerCmd = "docker run -w %s %s %s --rm %s %s"
      val dockerVolume = "-v %s:%s"
      val dockerVolumeInputs = s"$dockerVolume:ro"
      // `v.get` is safe below since we filtered the list earlier with only defined elements
      val inputVolumes = dockerInputDataVol.distinct.map(v => dockerVolumeInputs.format(v, v)).mkString(" ")
      val outputVolume = dockerVolume.format(executionDir.toAbsolutePath.toString, runtimeAttributes.dockerOutputDir.getOrElse(executionDir.toAbsolutePath.toString))
      val cmd = dockerCmd.format(runtimeAttributes.dockerWorkingDir.getOrElse(executionDir.toAbsolutePath.toString), inputVolumes, outputVolume, runtimeAttributes.dockerImage.get, jobCmd.get)
      log.debug("{} Docker command line to be used for task execution: {}.", tag, cmd)
      cmd
    }
  }

  private def prepareAndExecute: Unit = {
    Try {
      createExecutionFolderAndScript()
      executionResponse success executeTask()
    } recover {
      case exception => executionResponse success FailedNonRetryableResponse(jobDescriptor.key, exception, None)
    }
  }
}
