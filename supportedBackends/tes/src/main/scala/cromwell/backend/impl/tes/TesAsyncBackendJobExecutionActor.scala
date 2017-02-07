package cromwell.backend.impl.tes

import java.nio.file.FileAlreadyExistsException

import cromwell.backend.BackendJobLifecycleActor
import cromwell.backend.async.{ExecutionHandle, FailedNonRetryableExecutionHandle, PendingExecutionHandle}
import cromwell.backend.impl.tes.TesResponseJsonFormatter._
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.core.path.Obsolete._
import cromwell.core.path.Path
import cromwell.core.retry.SimpleExponentialBackoff
import spray.client.pipelining._
import spray.http.HttpRequest
import spray.httpx.SprayJsonSupport._
import spray.httpx.unmarshalling._
import wdl4s.values.WdlFile

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success}

sealed trait TesRunStatus {
  def isTerminal: Boolean
}

case object Running extends TesRunStatus {
  def isTerminal = false
}

case object Complete extends TesRunStatus {
  def isTerminal = true
}

case object FailedOrError extends TesRunStatus {
  def isTerminal = true
}

object TesAsyncBackendJobExecutionActor {
  val JobIdKey = "tes_job_id"
}

class TesAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
  extends BackendJobLifecycleActor with StandardAsyncExecutionActor with TesJobCachingActorHelper {

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

  override lazy val retryable: Boolean = false

  private val tesEndpoint = tesConfiguration.endpointURL

  override lazy val jobTag = jobDescriptor.key.tag

  private def pipeline[T: FromResponseUnmarshaller]: HttpRequest => Future[T] = sendReceive ~> unmarshal[T]

  // Utility for converting a WdlValue so that the path is localized to the
  // container's filesystem.
  override def mapCommandLineWdlFile(wdlFile: WdlFile): WdlFile = {
    val localPath = Paths.get(wdlFile.valueString).toAbsolutePath
    // Workaround since each input seemed to hit this function twice?
    if (localPath.startsWith(tesJobPaths.callInputsDockerRoot)) {
      WdlFile(localPath.pathAsString)
    } else {
      val containerPath = tesJobPaths.containerInput(localPath.pathAsString)
      WdlFile(containerPath)
    }
  }

  override lazy val commandDirectory: Path = tesJobPaths.callExecutionDockerRoot

  def createTaskMessage(): TesTaskMessage = {
    tesJobPaths.script.write(commandScriptContents)

    val task = TesTask(jobDescriptor, configurationDescriptor, jobLogger, tesJobPaths, runtimeAttributes)

    TesTaskMessage(
      Option(task.name),
      Option(task.description),
      Option(task.project),
      Option(task.inputs),
      Option(task.outputs),
      task.resources,
      task.dockerExecutor
    )
  }

  override def executeAsync()(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    // create call exec dir
    File(tesJobPaths.callExecutionRoot).createPermissionedDirectories()
    val taskMessage = createTaskMessage()

    val submitTask = pipeline[TesPostResponse]
      .apply(Post(tesEndpoint, taskMessage))

    submitTask.map {
      response =>
        val jobID = response.value
        PendingExecutionHandle(jobDescriptor, StandardAsyncJob(jobID), None, previousStatus = None)
    }
  }

  override def tryAbort(job: StandardAsyncJob): Unit = {

    val returnCodeTmp = jobPaths.returnCode.plusExt("kill")
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

  override def pollStatusAsync(handle: StandardAsyncPendingExecutionHandle)(implicit ec: ExecutionContext): Future[TesRunStatus] = {
    val pollTask = pipeline[TesGetResponse].apply(Get(s"$tesEndpoint/${handle.pendingJob.jobId}"))

    pollTask.map {
      response =>
        val state = response.state
        state match {
          case s if s.contains("Complete") => {
            jobLogger.info(s"Job ${handle.pendingJob.jobId} is complete")
            Complete
          }
          case s if s.contains("Cancel") => {
            jobLogger.info(s"Job ${handle.pendingJob.jobId} was canceled")
            FailedOrError
          }
          case s if s.contains("Error") => {
            jobLogger.info(s"TES reported an error for Job ${handle.pendingJob.jobId}")
            FailedOrError
          }
          case _ => Running
        }
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

  override def isSuccess(runStatus: TesRunStatus): Boolean = {
    runStatus match {
      case Complete => true
      case _ => false
    }
  }

  // Everything below was 'borrowed' from the SFS backend
  private def hostAbsoluteFilePath(wdlFile: WdlFile): File = {
    jobPaths.callRoot.resolve(wdlFile.value).toAbsolutePath
  }

  override def mapOutputWdlFile(wdlFile: WdlFile): WdlFile = {
    if (!hostAbsoluteFilePath(wdlFile).exists) {
      throw new RuntimeException("Could not process output, file not found: " +
        s"${hostAbsoluteFilePath(wdlFile).pathAsString}")
    } else {
      WdlFile(hostAbsoluteFilePath(wdlFile).pathAsString)
    }
  }
}
