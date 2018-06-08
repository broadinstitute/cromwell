package cromwell.backend.impl.tes

import java.io.FileNotFoundException
import java.nio.file.FileAlreadyExistsException

import akka.http.scaladsl.Http
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.marshalling.Marshal
import akka.http.scaladsl.model._
import akka.http.scaladsl.unmarshalling.{Unmarshal, Unmarshaller}
import akka.stream.ActorMaterializer
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.BackendJobLifecycleActor
import cromwell.backend.async.{ExecutionHandle, FailedNonRetryableExecutionHandle, PendingExecutionHandle}
import cromwell.backend.impl.tes.TesResponseJsonFormatter._
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.core.retry.SimpleExponentialBackoff
import wom.values.{WomFile, WomSingleFile}

import scala.concurrent.Future
import scala.concurrent.duration._
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
  implicit val actorSystem = context.system
  implicit val materializer = ActorMaterializer()

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

  // Utility for converting a WomValue so that the path is localized to the
  // container's filesystem.
  override def mapCommandLineWomFile(womFile: WomFile): WomFile = {
    womFile mapFile { path =>
      val localPath = DefaultPathBuilder.get(path).toAbsolutePath
      localPath match {
        case p if p.startsWith(tesJobPaths.workflowPaths.DockerRoot) => p.pathAsString
        case p if p.startsWith(tesJobPaths.callExecutionRoot) =>
          tesJobPaths.containerExec(commandDirectory, localPath.getFileName.pathAsString)
        case p if p.startsWith(tesJobPaths.callRoot) =>
          tesJobPaths.callDockerRoot.resolve(localPath.getFileName.pathAsString).pathAsString
        case p => tesJobPaths.containerInput(p.pathAsString)
      }
    }
  }

  override lazy val commandDirectory: Path = {
    runtimeAttributes.dockerWorkingDir match {
      case Some(path) => DefaultPathBuilder.get(path)
      case None => tesJobPaths.callExecutionDockerRoot
    }
  }

  def createTaskMessage(): ErrorOr[Task] = {
    val task =
      commandScriptContents.map(
        TesTask(
          jobDescriptor,
          configurationDescriptor,
          jobLogger,
          tesJobPaths,
          runtimeAttributes,
          commandDirectory,
          _,
          instantiatedCommand,
          realDockerImageUsed))

    task.map(task => Task(
      id = None,
      state = None,
      name = Option(task.name),
      description = Option(task.description),
      inputs = Option(task.inputs(commandLineValueMapper)),
      outputs = Option(task.outputs),
      resources = Option(task.resources),
      executors = task.executors,
      volumes = None,
      tags = None,
      logs = None
    ))
  }

  override def executeAsync(): Future[ExecutionHandle] = {
    // create call exec dir
    tesJobPaths.callExecutionRoot.createPermissionedDirectories()
    val taskMessageFuture = createTaskMessage().fold(
      errors => Future.failed(new RuntimeException(errors.toList.mkString(", "))),
      Future.successful)

    for {
      taskMessage <- taskMessageFuture
      entity <- Marshal(taskMessage).to[RequestEntity]
      ctr <- makeRequest[CreateTaskResponse](HttpRequest(method = HttpMethods.POST, uri = tesEndpoint, entity = entity))
    } yield PendingExecutionHandle(jobDescriptor, StandardAsyncJob(ctr.id), None, previousStatus = None)
  }

  override def recoverAsync(jobId: StandardAsyncJob) = executeAsync()

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

    makeRequest[CancelTaskResponse](HttpRequest(method = HttpMethods.POST, uri = s"$tesEndpoint/${job.jobId}:cancel")) onComplete {
      case Success(_) => jobLogger.info("{} Aborted {}", tag: Any, job.jobId)
      case Failure(ex) => jobLogger.warn("{} Failed to abort {}: {}", tag, job.jobId, ex.getMessage)
    }

    ()
  }

  override def pollStatusAsync(handle: StandardAsyncPendingExecutionHandle): Future[TesRunStatus] = {
    makeRequest[MinimalTaskView](HttpRequest(uri = s"$tesEndpoint/${handle.pendingJob.jobId}?view=MINIMAL")) map {
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
            jobLogger.info(s"TES reported an error for Job ${handle.pendingJob.jobId}: '$s'")
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

  override def isDone(runStatus: TesRunStatus): Boolean = {
    runStatus match {
      case Complete => true
      case _ => false
    }
  }
  
  private val outputWomFiles: Seq[WomSingleFile] = {
    Seq.empty
    // TODO WOM: fix
//    jobDescriptor.call.task
//      .findOutputFiles(jobDescriptor.fullyQualifiedInputs, NoFunctions)
//      .filter(o => !DefaultPathBuilder.get(o.valueString).isAbsolute)
  }

  override def mapOutputWomFile(womFile: WomFile): WomFile = {
    womFile mapFile { path =>
      val absPath: Path = tesJobPaths.callExecutionRoot.resolve(path)
      absPath match {
        case _ if !absPath.exists && outputWomFiles.map(_.value).contains(path) =>
          throw new FileNotFoundException("Could not process output, file not found: " +
            s"${absPath.pathAsString}")
        case _ => absPath.pathAsString
      }
    }
  }

  private def makeRequest[A](request: HttpRequest)(implicit um: Unmarshaller[ResponseEntity, A]): Future[A] = {
    for {
      response <- Http().singleRequest(request)
      data <- Unmarshal(response.entity).to[A]
    } yield data
  }
}
