package cromwell.backend.impl.bcs

import java.io.FileNotFoundException

import better.files.File.OpenOptions
import com.aliyuncs.batchcompute.main.v20151111.BatchComputeClient
import com.aliyuncs.exceptions.{ClientException, ServerException}
import cromwell.backend._
import cromwell.backend.async.{ExecutionHandle, FailedNonRetryableExecutionHandle, PendingExecutionHandle}
import cromwell.backend.impl.bcs.RunStatus.{Finished, TerminalRunStatus}
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.core.path.{DefaultPathBuilder, Path, PathFactory}
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.core.ExecutionEvent
import cromwell.filesystems.oss.OssPath
import wom.callable.Callable.OutputDefinition
import wom.core.FullyQualifiedName
import wom.expression.NoIoFunctionSet
import wom.types.WomSingleFileType
import wom.values._

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.util.{Success, Try}

object BcsAsyncBackendJobExecutionActor {
  val JobIdKey = "__bcs_job_id"
}

final class BcsAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
  extends BackendJobLifecycleActor with StandardAsyncExecutionActor with BcsJobCachingActorHelper {

  type BcsPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, BcsJob, RunStatus]

  override type StandardAsyncRunInfo = BcsJob

  override type StandardAsyncRunStatus = RunStatus

  override lazy val pollBackOff = SimpleExponentialBackoff(1.second, 5.minutes, 1.1)

  override lazy val executeOrRecoverBackOff = SimpleExponentialBackoff(3.seconds, 30.seconds, 1.1)

  // override lazy val dockerImageUsed: Option[String] = runtimeAttributes.docker map {docker => docker.image}
  override lazy val dockerImageUsed: Option[String] = None
  override lazy val commandDirectory: Path = BcsJobPaths.BcsCommandDirectory.resolve(bcsJobPaths.callExecutionRoot.pathWithoutScheme)

  private[bcs] lazy val userTag = runtimeAttributes.tag.getOrElse("cromwell")
  private[bcs] lazy val jobName: String =
    List(userTag, jobDescriptor.workflowDescriptor.id.shortString, jobDescriptor.taskCall.identifier.localName.value)
      .mkString("_")
      // Avoid "Name ... must only contain characters within [a-zA-Z0-9_-] and not start with [0-9]."
      .replaceAll("[^a-zA-Z0-9_-]", "_")

  override lazy val jobTag: String = jobDescriptor.key.tag

  private lazy val bcsWorkflowInputMount: BcsMount = bcsWorkflowPaths.getWorkflowInputMounts
  private lazy val userDefinedMounts: List[BcsMount] = runtimeAttributes.mounts.toList.flatten :+ bcsWorkflowInputMount
  // TODO: With a bit of refactoring this mutable var can be converted to a def or lazy val
  private var inputMounts: List[BcsMount] = List.empty

  private[bcs] def ossPathToMount(ossPath: OssPath): BcsInputMount = {
    val tmp = DefaultPathBuilder.get("/" + ossPath.pathWithoutScheme)
    val dir = tmp.getParent
    val local = BcsJobPaths.BcsTempInputDirectory.resolve(dir.pathAsString.md5SumShort).resolve(tmp.getFileName)
    val ret = BcsInputMount(ossPath, local, writeSupport = false)
    if (!inputMounts.exists(mount => mount.src == ossPath && mount.dest == local)) {
      inputMounts :+= ret
    }

    ret
  }

  private[bcs] def womFileToMount(file: WomFile): Option[BcsInputMount] = file match {
    case path if userDefinedMounts exists(bcsMount => path.valueString.startsWith(bcsMount.src.pathAsString)) => None
    case path => PathFactory.buildPath(path.valueString, initializationData.pathBuilders) match {
      case ossPath: OssPath => Some(ossPathToMount(ossPath))
      case _ => None
    }
  }

  private def bcsInputsFromWomFiles(prefix: String,
                                    remotePathArray: Seq[WomFile],
                                    jobDescriptor: BackendJobDescriptor): Iterable[BcsInputMount] = {
    remotePathArray flatMap { remotePath => womFileToMount(remotePath) match {
        case Some(mount) => Seq(mount)
        case None => Seq.empty
      }
    }
  }

  private[bcs] def getInputFiles(jobDescriptor: BackendJobDescriptor): Map[FullyQualifiedName, Seq[WomFile]] = {
    val writeFunctionFiles = instantiatedCommand.createdFiles map { f => f.file.value.md5SumShort -> Seq(f.file) }

    val writeFunctionInputs = writeFunctionFiles map  {
      case (name, files) => name -> files
    }

    // Collect all WomFiles from inputs to the call.
    val callInputFiles: Map[FullyQualifiedName, Seq[WomFile]] = jobDescriptor.fullyQualifiedInputs mapValues {
      _.collectAsSeq { case w: WomFile => w }
    }

    callInputFiles ++ writeFunctionInputs
  }

  private[bcs] def generateBcsInputs(jobDescriptor: BackendJobDescriptor): Unit = {
    val _ = getInputFiles(jobDescriptor) flatMap {
      case (name, files) => bcsInputsFromWomFiles(name, files, jobDescriptor)
    }
  }

  private def relativePath(path: String): Path = {
    val absolutePath = DefaultPathBuilder.get(path) match {
      case p if !p.isAbsolute => commandDirectory.resolve(p)
      case p => p
    }

    absolutePath
  }

  private[bcs] lazy val callRawOutputFiles: List[WomFile] = {
    import cats.syntax.validated._
    def evaluateFiles(output: OutputDefinition): List[WomFile] = {
      Try (
        output.expression.evaluateFiles(jobDescriptor.localInputs, NoIoFunctionSet, output.womType).map(_.toList map { _.file })
      ).getOrElse(List.empty[WomFile].validNel)
        .getOrElse(List.empty)
    }

    // val womFileOutputs = call.task.findOutputFiles(jobDescriptor.fullyQualifiedInputs, PureStandardLibraryFunctions)

    jobDescriptor.taskCall.callable.outputs.flatMap(evaluateFiles)
  }

  private[bcs] def isOutputOssFileString(s: String): Boolean = {
    callRawOutputFiles.exists({
      case file: WomSingleFile if file.value == s => true
      case _ => false
    })
  }

  private[bcs] def generateBcsOutputs(jobDescriptor: BackendJobDescriptor): Seq[BcsMount] = {
    callRawOutputFiles.flatMap(_.flattenFiles).distinct flatMap { womFile =>
      womFile match {
        case singleFile: WomSingleFile => List(generateBcsSingleFileOutput(singleFile))
        case _: WomGlobFile => throw new RuntimeException(s"glob output not supported currently")
        case _: WomUnlistedDirectory => throw new RuntimeException(s"directory output not supported currently")
      }
    }
  }

  private def generateBcsSingleFileOutput(wdlFile: WomSingleFile): BcsOutputMount = {
    val destination = getPath(wdlFile.valueString) match {
      case Success(ossPath: OssPath) => ossPath
      case Success(path: Path) if !path.isAbsolute => relativeOutputPath(path)
      case _ => callRootPath.resolve(wdlFile.value.stripPrefix("/"))
    }

    val src = relativePath(wdlFile.valueString)

    BcsOutputMount(src, destination, writeSupport = false)
  }

  private[bcs] def getOssFileName(ossPath: OssPath): String = {
    getPath(ossPath.pathWithoutScheme) match {
      case Success(path) => path.getFileName.pathAsString
      case _ => ossPath.pathWithoutScheme
    }
  }

  private[bcs] def localizeOssPath(ossPath: OssPath): String = {
    if (isOutputOssFileString(ossPath.pathAsString) && !ossPath.isAbsolute) {
      if (ossPath.exists) {
        ossPathToMount(ossPath).dest.normalize().pathAsString
      } else {
        commandDirectory.resolve(getOssFileName(ossPath)).normalize().pathAsString
      }
    } else {
      userDefinedMounts collectFirst {
        case bcsMount: BcsMount if ossPath.pathAsString.startsWith(bcsMount.src.pathAsString) =>
          bcsMount.dest.resolve(ossPath.pathAsString.stripPrefix(bcsMount.src.pathAsString)).pathAsString
      } getOrElse {
        val mount = ossPathToMount(ossPath)
        mount.dest.pathAsString
      }
    }
  }

  private[bcs] def relativeOutputPath(path: Path): Path = {
    if (isOutputOssFileString(path.pathAsString)) {
      bcsJobPaths.callRoot.resolve(path.pathAsString).normalize()
    } else {
      path
    }
  }

  override def mapCommandLineWomFile(womFile: WomFile): WomFile = {
    getPath(womFile.valueString) match {
      case Success(ossPath: OssPath) =>
        WomFile(WomSingleFileType, localizeOssPath(ossPath))
      case Success(path: Path) if !path.isAbsolute =>
        WomFile(WomSingleFileType, relativeOutputPath(path).pathAsString)
      case _ => womFile
    }
  }

  override def isTerminal(runStatus: RunStatus): Boolean = {
    runStatus match {
      case _ : TerminalRunStatus => true
      case _ => false
    }
  }

  override def getTerminalEvents(runStatus: RunStatus): Seq[ExecutionEvent] = {
    runStatus match {
      case successStatus: Finished => successStatus.eventList
      case unknown =>
        throw new RuntimeException(s"handleExecutionSuccess not called with RunStatus.Success. Instead got $unknown")
    }
  }

  override def handleExecutionFailure(runStatus: RunStatus,
                                      returnCode: Option[Int]): Future[ExecutionHandle] = {
    runStatus match {
      case RunStatus.Failed(jobId, Some(errorMessage), _) =>
        val exception = new Exception(s"Job id $jobId failed: '$errorMessage'")
        Future.successful(FailedNonRetryableExecutionHandle(exception, returnCode))
      case _ => super.handleExecutionFailure(runStatus, returnCode)
    }
  }

  override def isDone(runStatus: RunStatus): Boolean = {
    runStatus match {
      case _: Finished =>
        runtimeAttributes.autoReleaseJob match {
          case Some(true) | None =>
            bcsClient.deleteJob(runStatus.jobId)
          case _ =>
        }
        true
      case _ => false
    }
  }

  private[bcs] lazy val rcBcsOutput = BcsOutputMount(
    commandDirectory.resolve(bcsJobPaths.returnCodeFilename), bcsJobPaths.returnCode, writeSupport = false)

  private[bcs] lazy val stdoutBcsOutput = BcsOutputMount(
    commandDirectory.resolve(bcsJobPaths.defaultStdoutFilename), standardPaths.output, writeSupport = false)
  private[bcs] lazy val stderrBcsOutput = BcsOutputMount(
    commandDirectory.resolve(bcsJobPaths.defaultStderrFilename), standardPaths.error, writeSupport = false)

  private[bcs] lazy val uploadBcsWorkerPackage = {
    getPath(runtimeAttributes.workerPath.getOrElse(bcsJobPaths.workerFileName)) match {
      case Success(ossPath: OssPath) =>
        if (ossPath.notExists) {
          throw new FileNotFoundException(s"$ossPath")
        }
        ossPath
      case Success(path: Path) =>
        if (path.notExists) {
          throw new FileNotFoundException(s"$path")
        }

        if (bcsJobPaths.workerPath.notExists) {
          val content = path.byteArray
          bcsJobPaths.workerPath.writeByteArray(content)(OpenOptions.default)
        }

        bcsJobPaths.workerPath
      case _ => throw new RuntimeException(s"Invalid worker packer path")
    }
  }

  override def executeAsync(): Future[ExecutionHandle] = {
    commandScriptContents.fold(
      errors => Future.failed(new RuntimeException(errors.toList.mkString(", "))),
      bcsJobPaths.script.write)


    setBcsVerbose()

    val bcsJob = new BcsJob(
          jobName,
          jobTag,
          bcsCommandLine,
          uploadBcsWorkerPackage,
          bcsMounts,
          bcsEnvs,
          runtimeAttributes,
          Some(bcsJobPaths.bcsStdoutPath),
          Some(bcsJobPaths.bcsStderrPath),
          bcsClient)

    for {
      jobId <- Future.fromTry(bcsJob.submit())
    } yield PendingExecutionHandle(jobDescriptor, StandardAsyncJob(jobId), Option(bcsJob), previousStatus = None)
  }

  override def recoverAsync(jobId: StandardAsyncJob) = executeAsync()

  override def pollStatusAsync(handle: BcsPendingExecutionHandle): Future[RunStatus] = {
    val jobId = handle.pendingJob.jobId
    val bcsJob: BcsJob = handle.runInfo.getOrElse(throw new RuntimeException("empty run job info "))

    Future.fromTry(bcsJob.getStatus(jobId))
  }

  override def mapOutputWomFile(wdlFile: WomFile): WomFile = {
    wdlFileToOssPath(generateBcsOutputs(jobDescriptor))(wdlFile)
  }

  private[bcs] def wdlFileToOssPath(bcsOutputs: Seq[BcsMount])(wdlFile: WomFile): WomFile = {
    bcsOutputs collectFirst {
      case bcsOutput if bcsOutput.src.pathAsString.endsWith(wdlFile.valueString) => WomFile(WomSingleFileType, bcsOutput.dest.pathAsString)
    } getOrElse wdlFile
  }

  override def tryAbort(job: StandardAsyncJob): Unit = {
    for {
      client <- Try(initializationData.bcsConfiguration.bcsClient getOrElse(throw new RuntimeException("empty run job info ")))
      resp <- Try(client.getJob(job.jobId))
      status <- RunStatusFactory.getStatus(job.jobId, resp.getJob.getState)
    } yield {
      status match {
        case _: RunStatus.TerminalRunStatus =>
          for {
            _ <- Try(client.deleteJob(job.jobId))
          } yield job
        case _ =>
          for {
            _ <- Try(client.stopJob(job.jobId))
            _ <- Try(client.deleteJob(job.jobId))
          } yield job
      }
    }
    ()
  }

  override def isFatal(throwable: Throwable): Boolean = super.isFatal(throwable) || isFatalBcsException(throwable)

  private[bcs] def isFatalBcsException(throwable: Throwable) = {
    throwable match {
      case e: ClientException if e.getErrCode.startsWith("Invalid") => true
      case _ => false
    }
  }

  override def isTransient(throwable: Throwable): Boolean = {
    throwable match {
      case _: ServerException => true
      case e: ClientException if e.getErrCode == "InternalError" => true
      case _ => false
    }
  }

  private[bcs] def setBcsVerbose(): Unit = {
    runtimeAttributes.verbose match {
      case Some(verbose) => BatchComputeClient.verbose = verbose
      case None => BatchComputeClient.verbose = false
    }
  }

  private[bcs] lazy val bcsEnvs: Map[String, String] = Map(
      BcsJobPaths.BcsEnvCwdKey -> commandDirectory.pathAsString,
      BcsJobPaths.BcsEnvExecKey -> bcsJobPaths.script.pathAsString,
      BcsJobPaths.BcsEnvStdoutKey -> commandDirectory.resolve(bcsJobPaths.defaultStdoutFilename).pathAsString,
      BcsJobPaths.BcsEnvStderrKey -> commandDirectory.resolve(bcsJobPaths.defaultStderrFilename).pathAsString,
      BcsConfiguration.OssEndpointKey -> bcsConfiguration.ossEndpoint,
      BcsConfiguration.OssIdKey -> bcsConfiguration.ossAccessId,
      BcsConfiguration.OssSecretKey -> bcsConfiguration.ossAccessKey,
      BcsConfiguration.OssTokenKey -> bcsConfiguration.ossSecurityToken
  )

  private[bcs] lazy val bcsMounts: Seq[BcsMount] ={
    generateBcsInputs(jobDescriptor)
    runtimeAttributes.mounts.getOrElse(Seq.empty) ++ inputMounts ++
      generateBcsOutputs(jobDescriptor) :+ rcBcsOutput :+ stdoutBcsOutput :+ stderrBcsOutput :+ bcsWorkflowInputMount
  }
}
