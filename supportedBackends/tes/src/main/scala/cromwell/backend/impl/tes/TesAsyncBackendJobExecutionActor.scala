package cromwell.backend.impl.tes

import java.nio.file.FileAlreadyExistsException

import cromwell.backend.BackendJobLifecycleActor
import cromwell.backend.async.{ExecutionHandle, FailedNonRetryableExecutionHandle, PendingExecutionHandle}
import cromwell.backend.impl.tes.TesResponseJsonFormatter._
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.core.retry.SimpleExponentialBackoff
import spray.client.pipelining._
import spray.http.HttpRequest
import spray.httpx.SprayJsonSupport._
import spray.httpx.unmarshalling._
import wdl4s.expression.NoFunctions
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
  
  private lazy val realDockerImageUsed: String = jobDescriptor.maybeCallCachingEligible.dockerHash.getOrElse(runtimeAttributes.dockerImage)
  override lazy val dockerImageUsed: Option[String] = Option(realDockerImageUsed)

  private val tesEndpoint = workflowDescriptor.workflowOptions.getOrElse("endpoint", tesConfiguration.endpointURL)

  override lazy val jobTag: String = jobDescriptor.key.tag

  private def pipeline[T: FromResponseUnmarshaller]: HttpRequest => Future[T] = sendReceive ~> unmarshal[T]

  // Utility for converting a WdlValue so that the path is localized to the
  // container's filesystem.
  override def mapCommandLineWdlFile(wdlFile: WdlFile): WdlFile = {
    val localPath = DefaultPathBuilder.get(wdlFile.valueString).toAbsolutePath
    localPath match {
      case p if p.startsWith(tesJobPaths.workflowPaths.DockerRoot) =>
        val containerPath = p.pathAsString
        WdlFile(containerPath)
      case p if p.startsWith(tesJobPaths.callExecutionRoot) =>
        val containerPath = tesJobPaths.containerExec(commandDirectory, localPath.getFileName.pathAsString)
        WdlFile(containerPath)
      case p =>
        val containerPath = tesJobPaths.containerInput(p.pathAsString)
        WdlFile(containerPath)
    }
  }

  override lazy val commandDirectory: Path = {
    runtimeAttributes.dockerWorkingDir match {
      case Some(path) => DefaultPathBuilder.get(path)
      case None => tesJobPaths.callExecutionDockerRoot
    }
  }

  def createTaskMessage(): Task = {
    val task = TesTask(jobDescriptor, configurationDescriptor, jobLogger, tesJobPaths,
      runtimeAttributes, commandDirectory, commandScriptContents, backendEngineFunctions,
      realDockerImageUsed)

    Task(
      None,
      None,
      Option(task.name),
      Option(task.description),
      Option(task.project),
      Option(task.inputs(commandLineValueMapper)),
      Option(task.outputs),
      Option(task.resources),
      task.executors,
      None,
      None,
      None
    )
  }

  override def executeAsync()(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    // create call exec dir
    tesJobPaths.callExecutionRoot.createPermissionedDirectories()
    val taskMessage = createTaskMessage()

    val submitTask = pipeline[CreateTaskResponse]
      .apply(Post(tesEndpoint, taskMessage))

    submitTask.map {
      response =>
        val jobID = response.id
        PendingExecutionHandle(jobDescriptor, StandardAsyncJob(jobID), None, previousStatus = None)
    }
  }

  override def recoverAsync(jobId: StandardAsyncJob)(implicit ec: ExecutionContext) = executeAsync()

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

    val abortRequest = pipeline[CancelTaskResponse]
      .apply(Post(s"$tesEndpoint/${job.jobId}:cancel"))
    abortRequest onComplete {
      case Success(_) => jobLogger.info("{} Aborted {}", tag: Any, job.jobId)
      case Failure(ex) => jobLogger.warn("{} Failed to abort {}: {}", tag, job.jobId, ex.getMessage)
    }
    ()
  }

  override def pollStatusAsync(handle: StandardAsyncPendingExecutionHandle)(implicit ec: ExecutionContext): Future[TesRunStatus] = {
    val pollTask = pipeline[MinimalTaskView].apply(Get(s"$tesEndpoint/${handle.pendingJob.jobId}?view=MINIMAL"))

    pollTask.map {
      response =>
        val state = response.state
        state match {
          case s if s.contains("COMPLETE") =>
            jobLogger.info(s"Job ${handle.pendingJob.jobId} is complete")
            Complete

          case s if s.contains("CANCELED") =>
            jobLogger.info(s"Job ${handle.pendingJob.jobId} was canceled")
            FailedOrError

          case s if s.contains("ERROR") =>
            jobLogger.info(s"TES reported an error for Job ${handle.pendingJob.jobId}")
            FailedOrError

          case _ => Running
        }
    }
  }

  override def customPollStatusFailure: PartialFunction[(ExecutionHandle, Exception), ExecutionHandle] = {
    case (oldHandle: StandardAsyncPendingExecutionHandle@unchecked, e: Exception) =>
      jobLogger.error(s"$tag TES Job ${oldHandle.pendingJob.jobId} has not been found, failing call")
      FailedNonRetryableExecutionHandle(e)
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
  
  private val outputWdlFiles: Seq[WdlFile] = jobDescriptor.call.task
    .findOutputFiles(jobDescriptor.fullyQualifiedInputs, NoFunctions)
    .filter(o => !DefaultPathBuilder.get(o.valueString).isAbsolute)

  override def mapOutputWdlFile(wdlFile: WdlFile): WdlFile = {
    val absPath: Path = tesJobPaths.callExecutionRoot.resolve(wdlFile.valueString)
    wdlFile match {
      case fileNotFound if !absPath.exists && outputWdlFiles.contains(fileNotFound) =>
        throw new RuntimeException("Could not process output, file not found: " +
          s"${absPath.pathAsString}")
      case _ => WdlFile(absPath.pathAsString)
    }
  }
}
