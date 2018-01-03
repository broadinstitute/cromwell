package cromwell.backend.impl.bcs



import java.io.FileNotFoundException

import akka.stream.ActorMaterializer
import com.aliyuncs.batchcompute.main.v20151111.BatchComputeClient
import com.aliyuncs.exceptions.{ClientException, ServerException}
import cromwell.backend._
import cromwell.backend.async.{ExecutionHandle, PendingExecutionHandle}
import cromwell.backend.impl.bcs.RunStatus.{Finished, TerminalRunStatus}
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.core.{ExecutionEvent, NoIoFunctionSet}
import cromwell.filesystems.oss.{OssPath, UploadFileOption, UploadStringOption}
import cromwell.core.retry.SimpleExponentialBackoff

import scala.concurrent.duration._
import cromwell.core.path.{DefaultPathBuilder, Path, PathFactory}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import wom.WomFileMapper
import wom.callable.Callable.OutputDefinition
import wom.core.FullyQualifiedName
import wom.types.WomSingleFileType
import wom.values._

import scala.concurrent.Future
import scala.util.{Success, Try}


object BcsAsyncBackendJobExecutionActor {
  val JobIdKey = "__bcs_job_id"
}

class BcsAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
  extends BackendJobLifecycleActor with StandardAsyncExecutionActor with BcsJobCachingActorHelper with GcsBatchCommandBuilder {

  implicit val actorSystem = context.system
  implicit val materializer = ActorMaterializer()


  type BcsPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, BcsJob, RunStatus]

  override type StandardAsyncRunInfo = BcsJob

  override type StandardAsyncRunStatus = RunStatus

  override lazy val pollBackOff = SimpleExponentialBackoff(1.second, 5.minutes, 1.1)

  override lazy val executeOrRecoverBackOff = SimpleExponentialBackoff(3.seconds, 30.seconds, 1.1)

  override lazy val dockerImageUsed: Option[String] = runtimeAttributes.docker map {docker => docker.image}
  override lazy val commandDirectory: Path = BcsJobPaths.BcsCommandDirectory.resolve(bcsJobPaths.callExecutionRoot.pathWithoutScheme)

  private[bcs] lazy val jobName: String = s"cromwell_${jobDescriptor.workflowDescriptor.id.shortString}_${jobDescriptor.taskCall.identifier.localName}"
  override lazy val jobTag: String = jobDescriptor.key.tag

  private lazy val bcsWorkflowInputMount: BcsMount = bcsWorkflowPaths.getWorkflowInputMounts
  private lazy val userDefinedMounts = runtimeAttributes.mounts.getOrElse(Seq.empty) :+ bcsWorkflowInputMount
  private[bcs] def localizeOssSingleInput(ossPath: OssPath): Path = {
    val tmp = DefaultPathBuilder.get("/" + ossPath.pathWithoutScheme)
    val dir = tmp.getParent
    BcsJobPaths.BcsTempInputDirectory.resolve(dir.pathAsString.md5SumShort).resolve(tmp.getFileName)
  }

  private def bcsInputsFromWomFiles(prefix: String,
                                    remotePathArray: Seq[WomFile],
                                    jobDescriptor: BackendJobDescriptor): Iterable[BcsInputMount] = {
    remotePathArray flatMap {
      case remotePath if userDefinedMounts exists(bcsMount => remotePath.valueString.startsWith(bcsMount.src.pathAsString)) => Seq.empty
      case remotePath => PathFactory.buildPath(remotePath.valueString, initializationData.pathBuilders) match {
        case ossPath: OssPath => Seq(BcsInputMount(ossPath, localizeOssSingleInput(ossPath), writeSupport = false))
        case _ => Seq.empty
      }
    }
  }

  private[bcs] def getInputFiles(jobDescriptor: BackendJobDescriptor): Map[FullyQualifiedName, Seq[WomFile]] = {
    /*
    val fullyQualifiedPreprocessedInputs = jobDescriptor.inputDeclarations map { case (declaration, value) => declaration.fullyQualifiedName -> value }
    val writeFunctionFiles = call.task.evaluateFilesFromCommand(fullyQualifiedPreprocessedInputs, backendEngineFunctions) map {
      case (expression, file) => expression.toWdlString.md5SumShort -> Seq(file)
    }

    val callInputFiles: Map[FullyQualifiedName, Seq[WomFile]] = jobDescriptor.fullyQualifiedInputs mapValues {
      _.collectAsSeq { case w: WomFile => w }
    }
    */
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

  private[bcs] lazy val inputFiles = getInputFiles(jobDescriptor)

  private[bcs] def isInputFile(pathAsString: String): Boolean = {
    inputFiles.values exists (
      files => files.exists(wdlFile => wdlFile.valueString == pathAsString)
    )
  }

  private[bcs] def generateBcsInputs(jobDescriptor: BackendJobDescriptor): Seq[BcsMount] = {
    val inputs = getInputFiles(jobDescriptor) flatMap {
      case (name, files) => bcsInputsFromWomFiles(name, files, jobDescriptor)
    }

    inputs.toSeq
  }

  private def relativePath(path: String): Path = {
    val absolutePath = DefaultPathBuilder.get(path) match {
      case p if !p.isAbsolute => commandDirectory.resolve(p)
      case p => p
    }

    absolutePath
  }

  private[bcs] def generateBcsOutputs(jobDescriptor: BackendJobDescriptor): Seq[BcsMount] = {
    import cats.syntax.validated._
    def evaluateFiles(output: OutputDefinition): List[WomFile] = {
      Try (
        output.expression.evaluateFiles(jobDescriptor.localInputs, NoIoFunctionSet, output.womType).map(_.toList)
      ).getOrElse(List.empty[WomFile].validNel)
        .getOrElse(List.empty)
    }

    // val womFileOutputs = call.task.findOutputFiles(jobDescriptor.fullyQualifiedInputs, PureStandardLibraryFunctions)

    val womFileOutputs = jobDescriptor.taskCall.callable.outputs.flatMap(evaluateFiles)

    val outputs = womFileOutputs.distinct flatMap { womFile =>
      womFile match {
        case singleFile: WomSingleFile => List(generateBcsSingleFileOutputs(singleFile))
        case _: WomGlobFile => throw new RuntimeException(s"glob output not supported currently")
      }
    }

    outputs.toSeq
  }

  private def generateBcsSingleFileOutputs(wdlFile: WomSingleFile): BcsOutputMount = {
    val destination = getPath(wdlFile.valueString) match {
      case Success(ossPath: OssPath) => ossPath
      case _ => callRootPath.resolve(wdlFile.value.stripPrefix("/"))
    }

    val src = relativePath(mapCommandLineWomFile(wdlFile).valueString)

    BcsOutputMount(src, destination, false)
  }

  private[bcs] def getOssFileName(ossPath: OssPath): String = {
    getPath(ossPath.pathWithoutScheme) match {
      case Success(path) => path.getFileName.pathAsString
      case _ => ossPath.pathWithoutScheme
    }
  }

  private[bcs] def localizeOssPath(ossPath: OssPath): String = {
    if (isInputFile(ossPath.pathAsString)) {
      userDefinedMounts collectFirst {
        case bcsMount: BcsMount if ossPath.pathAsString.startsWith(bcsMount.src.pathAsString) =>
          bcsMount.dest.resolve(ossPath.pathAsString.stripPrefix(bcsMount.src.pathAsString)).pathAsString
      } getOrElse {
        generateBcsInputs(jobDescriptor) collectFirst {
          case bcsMount if bcsMount.src.pathAsString == ossPath.pathAsString => bcsMount.dest.pathAsString
        } getOrElse commandDirectory.resolve(getOssFileName(ossPath)).pathAsString
      }
    } else {
      commandDirectory.resolve(getOssFileName(ossPath)).pathAsString
    }

  }

  override def mapCommandLineWomFile(womFile: WomFile): WomFile = {
    getPath(womFile.valueString) match {
      case Success(ossPath: OssPath) =>
        WomFile(WomSingleFileType, localizeOssPath(ossPath))
      case _ => womFile
    }
  }

  override lazy val commandLineValueMapper: WomValue => WomValue = {
    wdlValue => wdlValue match {
      case str: WomString => getPath(str.valueString) match {
        case Success(ossPath: OssPath) => WomString(localizeOssPath(ossPath))
        case _ => str
      }
      case wdlValue => WomFileMapper.mapWomFiles(mapCommandLineWomFile)(wdlValue).get
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
  override def isSuccess(runStatus: RunStatus): Boolean = {
    runStatus match {
      case _ : Finished => {
        runtimeAttributes.autoReleaseJob match {
          case Some(release) if release == true =>
            bcsClient.deleteJob(runStatus.jobId)
          case _ =>
        }
        true
      }
      case _ => false
    }
  }

  private[bcs] lazy val  rcBcsOutput = new BcsOutputMount(commandDirectory.resolve(bcsJobPaths.returnCodeFilename), bcsJobPaths.returnCode, false)
  private[bcs] lazy val  stdoutBcsOutput = new BcsOutputMount(commandDirectory.resolve(bcsJobPaths.stdoutFilename), bcsJobPaths.stdout, false)
  private[bcs] lazy val  stderrBcsOutput = new BcsOutputMount(commandDirectory.resolve(bcsJobPaths.stderrFilename), bcsJobPaths.stderr, false)

  private[bcs] lazy val uploadBcsWorkerPackage = {
    runtimeAttributes.workerPath match {
      case Some(pathAsString: String) =>
        getPath(pathAsString) match {
          case Success(ossPath: OssPath) =>
            if (ossPath.notExists) {
              throw new FileNotFoundException(s"$pathAsString")
            }
            ossPath
          case Success(path: Path) =>
            if (bcsJobPaths.workerPath.notExists) {
              writeAsync(bcsJobPaths.workerPath, path.pathAsString, Seq(UploadFileOption))
            }
            bcsJobPaths.workerPath
          case _ => throw new RuntimeException(s"Invalid worker packer path: $pathAsString")
        }
      case None =>
        writeAsync(bcsJobPaths.workerPath, bcsJobPaths.workerFileName, Seq(UploadFileOption))
        bcsJobPaths.workerPath
    }
  }

  override def executeAsync(): Future[ExecutionHandle] = {
    commandScriptContents.fold(
      errors => Future.failed(new RuntimeException(errors.toList.mkString(", "))),
      writeAsync(bcsJobPaths.script, _, Seq(UploadStringOption)))


    setBcsVerbose()

    val bcsJob = new BcsJob(
          jobName,
          jobTag,
          bcsCommandLine,
          uploadBcsWorkerPackage,
          bcsMounts,
          bcsEnvs,
          runtimeAttributes,
          bcsClient)

    for {
      jobId <- Future.fromTry(bcsJob.submit())
    } yield new PendingExecutionHandle(jobDescriptor, StandardAsyncJob(jobId), Option(bcsJob), previousStatus = None)
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
    return
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
      BcsJobPaths.BcsEnvStdoutKey -> commandDirectory.resolve(bcsJobPaths.stdoutFilename).pathAsString,
      BcsJobPaths.BcsEnvStderrKey -> commandDirectory.resolve(bcsJobPaths.stderrFilename).pathAsString,
      BcsConfiguration.OssEndpointKey -> bcsConfiguration.ossEndpoint,
      BcsConfiguration.OssIdKey -> bcsConfiguration.ossAccessId,
      BcsConfiguration.OssSecretKey -> bcsConfiguration.ossAccessKey,
      BcsConfiguration.OssTokenKey -> bcsConfiguration.ossSecurityToken
  )

  private[bcs] lazy val bcsMounts: Seq[BcsMount] =
    runtimeAttributes.mounts.getOrElse(Seq.empty) ++ generateBcsInputs(jobDescriptor) ++
      generateBcsOutputs(jobDescriptor) :+ rcBcsOutput :+ stdoutBcsOutput :+ stderrBcsOutput :+ bcsWorkflowInputMount

}
