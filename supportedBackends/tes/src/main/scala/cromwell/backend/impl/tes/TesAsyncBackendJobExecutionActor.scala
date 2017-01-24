package cromwell.backend.impl.tes

import java.nio.file.{FileAlreadyExistsException, Path}

import cromwell.backend.async.{ExecutionHandle, FailedNonRetryableExecutionHandle, PendingExecutionHandle, SuccessfulExecutionHandle}
import cromwell.backend.BackendJobLifecycleActor
import cromwell.core.retry.SimpleExponentialBackoff
import TesResponseJsonFormatter._
import better.files.File
import cromwell.backend.io.JobPaths
import cromwell.backend.sfs.SharedFileSystemExpressionFunctions
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob, StandardInitializationData}
import spray.httpx.SprayJsonSupport._
import spray.client.pipelining._
import spray.http.HttpRequest
import spray.httpx.unmarshalling._
import wdl4s.values.{WdlArray, WdlFile, WdlMap, WdlValue}
import cromwell.backend.wdl.OutputEvaluator
import cromwell.core.path.FileImplicits._
import cromwell.core.path.PathFactory.pathPlusSuffix
import cromwell.core.path.{PathBuilder, PathFactory}
import lenthall.util.TryUtil

import scala.concurrent.duration._
import scala.concurrent.{Await, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

case class TesRunStatus(isTerminal: Boolean)

object TesAsyncBackendJobExecutionActor {
  val JobIdKey = "tes_job_id"
}

class TesAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
  extends BackendJobLifecycleActor with StandardAsyncExecutionActor  with TesJobCachingActorHelper {

  override type StandardAsyncRunInfo = Any

  override type StandardAsyncRunStatus = TesRunStatus

  override lazy val pollBackOff = SimpleExponentialBackoff(
    initialInterval = 1 seconds,
    maxInterval = 5 minutes,
    multiplier = 1.1
  )

  override lazy val executeOrRecoverBackOff = SimpleExponentialBackoff(
    initialInterval = 3 seconds,
    maxInterval = 30 seconds,
    multiplier = 1.1
  )

  lazy val pathBuilders: List[PathBuilder] = StandardInitializationData.pathBuilders(backendInitializationDataOption)
  lazy val backendEngineFunctions = SharedFileSystemExpressionFunctions(jobPaths, pathBuilders)

  override lazy val retryable: Boolean = false

  private val tesEndpoint = tesConfiguration.endpointURL

  override lazy val jobTag = jobDescriptor.key.tag

  private def pipeline[T: FromResponseUnmarshaller]: HttpRequest => Future[T] = sendReceive ~> unmarshal[T]

  def createTaskMessage(): TesTaskMessage = {
    val task = TesTask(jobDescriptor, configurationDescriptor, jobLogger, tesJobPaths, runtimeAttributes)
    TesTaskMessage(
      Some(task.name),
      Some(task.description),
      Some(task.project),
      Some(task.inputs),
      Some(task.outputs),
      task.resources,
      task.dockerExecutor
    )
  }

  override def execute(): ExecutionHandle = {
    // create call exec dir
    File(tesJobPaths.callExecutionRoot).createPermissionedDirectories()
    val taskMessage = createTaskMessage()

    val submitTask = pipeline[TesPostResponse]
            .apply(Post(tesEndpoint, taskMessage))

    val response = Await.result(submitTask, 15.seconds)
    val jobID = response.value.get

    PendingExecutionHandle(jobDescriptor, StandardAsyncJob(jobID), None, previousStatus = None)
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

    val abortRequest = pipeline[TesGetResponse]
      .apply(Delete(s"$tesEndpoint/${job.jobId}"))
    abortRequest onComplete {
      case Success(_) => jobLogger.info("{} Aborted {}", tag: Any, job.jobId)
      case Failure(ex) => jobLogger.warn("{} Failed to abort {}: {}", tag, job.jobId, ex.getMessage)
    }
    ()
  }

  override def requestsAbortAndDiesImmediately: Boolean = true

  override def pollStatus(handle: StandardAsyncPendingExecutionHandle): TesRunStatus = {
    val pollTask = pipeline[TesGetResponse].apply(Get(s"$tesEndpoint/${handle.pendingJob.jobId}"))

    val response = Await.result(pollTask, 5.seconds)
    val state = response.state.get
    if (state contains "Complete") {
      jobLogger.info(s"Job ${handle.pendingJob.jobId} is complete")
      TesRunStatus(isTerminal = true)
    } else {
      TesRunStatus(isTerminal = false)
    }
  }

  override def customPollStatusFailure: PartialFunction[(ExecutionHandle, Exception), ExecutionHandle] = {
    case (oldHandle: StandardAsyncPendingExecutionHandle@unchecked, e: Exception) => {
      jobLogger.error(s"$tag TES Job ${oldHandle.pendingJob.jobId} has not been found, failing call")
      FailedNonRetryableExecutionHandle(e)
    }
  }

  override def isTerminal(runStatus: TesRunStatus): Boolean = {
    runStatus.isTerminal
  }

  // Everything below was 'borrowed' from the SFS backend
  private def hostAbsoluteFilePath(callRoot: Path, pathString: String): File = {
    val wdlPath = PathFactory.buildPath(pathString, pathBuilders)
    callRoot.resolve(wdlPath).toAbsolutePath
  }

  def outputMapper(job: JobPaths)(wdlValue: WdlValue): Try[WdlValue] = {
    wdlValue match {
      case fileNotFound: WdlFile if !hostAbsoluteFilePath(job.callExecutionRoot, fileNotFound.valueString).exists =>
        Failure(new RuntimeException("Could not process output, file not found: " +
          s"${hostAbsoluteFilePath(job.callExecutionRoot, fileNotFound.valueString).pathAsString}"))
      case file: WdlFile => Try(WdlFile(hostAbsoluteFilePath(job.callExecutionRoot, file.valueString).pathAsString))
      case array: WdlArray =>
        val mappedArray = array.value map outputMapper(job)
        TryUtil.sequence(mappedArray) map { WdlArray(array.wdlType, _) }
      case map: WdlMap =>
        val mappedMap = map.value mapValues outputMapper(job)
        TryUtil.sequenceMap(mappedMap) map { WdlMap(map.wdlType, _) }
      case other => Success(other)
    }
  }

  override def handleExecutionSuccess(runStatus: StandardAsyncRunStatus,
                                      handle: StandardAsyncPendingExecutionHandle,
                                      returnCode: Int): ExecutionHandle = {

    val outputsTry =
      OutputEvaluator.evaluateOutputs(jobDescriptor, backendEngineFunctions, outputMapper(jobPaths))
    outputsTry match {
      case Success(outputs) => {
        SuccessfulExecutionHandle(outputs, returnCode, jobPaths.detritusPaths, Seq.empty)
      }
      case Failure(throwable) => FailedNonRetryableExecutionHandle(throwable, Option(returnCode))
    }
  }
}
