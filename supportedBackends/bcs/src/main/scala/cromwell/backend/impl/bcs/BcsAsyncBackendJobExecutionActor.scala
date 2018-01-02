package cromwell.backend.impl.bcs



import java.io.FileNotFoundException

import akka.stream.ActorMaterializer
import com.aliyuncs.batchcompute.main.v20151111.BatchComputeClient
import com.aliyuncs.exceptions.{ClientException, ServerException}
import cromwell.backend._
import cromwell.backend.async.{ExecutionHandle, PendingExecutionHandle}
import cromwell.backend.impl.bcs.RunStatus.{Finished, TerminalRunStatus}
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.backend.wdl.WdlFileMapper
import cromwell.core.ExecutionEvent
import cromwell.filesystems.oss.{OssPath, UploadFileOption, UploadStringOption}
import cromwell.core.retry.SimpleExponentialBackoff
import wdl4s.wdl.values.WdlFile

import scala.concurrent.duration._
import cromwell.core.path.{DefaultPathBuilder, Path, PathFactory}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import wdl4s.wdl._
import wdl4s.wdl.expression.PureStandardLibraryFunctions
import wdl4s.wdl.values._

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

  private[bcs] lazy val jobName: String = s"cromwell_${jobDescriptor.workflowDescriptor.id.shortString}_${jobDescriptor.call.unqualifiedName}"
  override lazy val jobTag: String = jobDescriptor.key.tag

  private lazy val bcsWorkflowInputMount: BcsMount = bcsWorkflowPaths.getWorkflowInputMounts
  private lazy val userDefinedMounts = runtimeAttributes.mounts.getOrElse(Seq.empty) :+ bcsWorkflowInputMount
  private[bcs] def localizeOssSingleInput(ossPath: OssPath): Path = {
    val tmp = DefaultPathBuilder.get("/" + ossPath.pathWithoutScheme)
    val dir = tmp.getParent
    BcsJobPaths.BcsTempInputDirectory.resolve(dir.pathAsString.md5SumShort).resolve(tmp.getFileName)
  }

  private def bcsInputsFromWdlFiles(prefix: String,
                                    remotePathArray: Seq[WdlFile],
                                    jobDescriptor: BackendJobDescriptor): Iterable[BcsInputMount] = {
    remotePathArray flatMap {
      case remotePath if userDefinedMounts exists(bcsMount => remotePath.valueString.startsWith(bcsMount.src.pathAsString)) => Seq.empty
      case remotePath => PathFactory.buildPath(remotePath.valueString, initializationData.pathBuilders) match {
        case ossPath: OssPath => Seq(BcsInputMount(ossPath, localizeOssSingleInput(ossPath), writeSupport = false))
        case _ => Seq.empty
      }
    }
  }

  private[bcs] def getInputFiles(jobDescriptor: BackendJobDescriptor): Map[FullyQualifiedName, Seq[WdlFile]] = {
    val fullyQualifiedPreprocessedInputs = jobDescriptor.inputDeclarations map { case (declaration, value) => declaration.fullyQualifiedName -> value }
    val writeFunctionFiles = call.task.evaluateFilesFromCommand(fullyQualifiedPreprocessedInputs, backendEngineFunctions) map {
      case (expression, file) => expression.toWdlString.md5SumShort -> Seq(file)
    }

    val callInputFiles: Map[FullyQualifiedName, Seq[WdlFile]] = jobDescriptor.fullyQualifiedInputs mapValues {
      _.collectAsSeq { case w: WdlFile => w }
    }

    callInputFiles ++ writeFunctionFiles
  }

  private[bcs] lazy val inputFiles = getInputFiles(jobDescriptor)

  private[bcs] def isInputFile(pathAsString: String): Boolean = {
    inputFiles.values exists (
      files => files.exists(wdlFile => wdlFile.valueString == pathAsString)
    )
  }

  private[bcs] def generateBcsInputs(jobDescriptor: BackendJobDescriptor): Seq[BcsMount] = {
    val inputs = getInputFiles(jobDescriptor) flatMap {
      case (name, files) => bcsInputsFromWdlFiles(name, files, jobDescriptor)
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
    val wdlFileOutputs = call.task.findOutputFiles(jobDescriptor.fullyQualifiedInputs, PureStandardLibraryFunctions)

    val outputs = wdlFileOutputs.distinct flatMap { wdlFile =>
      wdlFile match {
        case singleFile: WdlSingleFile => List(generateBcsSingleFileOutputs(singleFile))
        case _: WdlGlobFile => throw new RuntimeException(s"glob output not supported currently")
      }
    }

    outputs.toSeq
  }

  private def generateBcsSingleFileOutputs(wdlFile: WdlSingleFile): BcsOutputMount = {
    val destination = getPath(wdlFile.valueString) match {
      case Success(ossPath: OssPath) => ossPath
      case _ => callRootPath.resolve(wdlFile.value.stripPrefix("/"))
    }

    val src = relativePath(mapCommandLineWdlFile(wdlFile).valueString)

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

  override def mapCommandLineWdlFile(wdlFile: WdlFile): WdlFile = {
    getPath(wdlFile.valueString) match {
      case Success(ossPath: OssPath) =>
        WdlFile(localizeOssPath(ossPath), wdlFile.isGlob)
      case _ => wdlFile
    }
  }

  override lazy val commandLineValueMapper: WdlValue => WdlValue = {
    wdlValue => wdlValue match {
      case str: WdlString => getPath(str.valueString) match {
        case Success(ossPath: OssPath) => WdlString(localizeOssPath(ossPath))
        case _ => str
      }
      case wdlValue => WdlFileMapper.mapWdlFiles(mapCommandLineWdlFile)(wdlValue).get
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
    writeAsync(bcsJobPaths.script, commandScriptContents, Seq(UploadStringOption))

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

  override def mapOutputWdlFile(wdlFile: WdlFile): WdlFile = {
    wdlFileToOssPath(generateBcsOutputs(jobDescriptor))(wdlFile)
  }

  private[bcs] def wdlFileToOssPath(bcsOutputs: Seq[BcsMount])(wdlFile: WdlFile): WdlFile = {
    bcsOutputs collectFirst {
      case bcsOutput if bcsOutput.src.pathAsString.endsWith(wdlFile.valueString) => WdlFile(bcsOutput.dest.pathAsString)
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
