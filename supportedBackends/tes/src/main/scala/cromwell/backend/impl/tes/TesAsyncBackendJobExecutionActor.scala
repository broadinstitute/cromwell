package cromwell.backend.impl.tes

import akka.event.LoggingAdapter
import akka.http.scaladsl.Http
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.marshalling.Marshal
import akka.http.scaladsl.model.HttpHeader.ParsingResult.Ok
import akka.http.scaladsl.model._
import akka.http.scaladsl.unmarshalling.{Unmarshal, Unmarshaller}
import akka.stream.ActorMaterializer
import akka.util.ByteString
import cats.implicits._
import common.exception.AggregatedMessageException
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cromwell.backend.async.{
  AbortedExecutionHandle,
  ExecutionHandle,
  FailedNonRetryableExecutionHandle,
  PendingExecutionHandle
}
import cromwell.backend.impl.tes.TesAsyncBackendJobExecutionActor._
import cromwell.backend.impl.tes.TesResponseJsonFormatter._
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.backend.{BackendJobLifecycleActor, Platform}
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.core.retry.Retry._
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.filesystems.drs.{DrsPath, DrsResolver}
import cromwell.filesystems.http.HttpPath
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.metadata.CallMetadataKeys
import net.ceedubs.ficus.Ficus._
import wom.values.WomFile

import java.io.FileNotFoundException
import java.nio.file.FileAlreadyExistsException
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

case class TesVmCostData(startTime: Option[String], endTime: Option[String], vmCost: Option[String]) {
  val fullyPopulated: Boolean = startTime.nonEmpty && vmCost.nonEmpty
}
sealed trait TesRunStatus {
  def isTerminal: Boolean
  def sysLogs: Seq[String] = Seq.empty[String]
  def costData: Option[TesVmCostData] = None
}

case class Running(override val costData: Option[TesVmCostData] = Option.empty) extends TesRunStatus {
  def isTerminal = false
  override def toString = "Running"
}

case class Complete(override val costData: Option[TesVmCostData] = Option.empty) extends TesRunStatus {
  def isTerminal = true
  override def toString = "Complete"
}

case class Error(override val sysLogs: Seq[String] = Seq.empty[String],
                 override val costData: Option[TesVmCostData] = Option.empty
) extends TesRunStatus {
  def isTerminal = true
  override def toString = "SYSTEM_ERROR"
}

case class Failed(override val sysLogs: Seq[String] = Seq.empty[String],
                  override val costData: Option[TesVmCostData] = Option.empty
) extends TesRunStatus {
  def isTerminal = true
  override def toString = "EXECUTOR_ERROR"
}

case class Cancelled(override val costData: Option[TesVmCostData] = Option.empty) extends TesRunStatus {
  def isTerminal = true
  override def toString = "Cancelled"
}

object TesAsyncBackendJobExecutionActor {
  val JobIdKey = "tes_job_id"

  type StandardAsyncRunInfo = Any
  type StandardAsyncRunState = TesRunStatus

  type StandardAsyncPendingExecutionHandle =
    PendingExecutionHandle[StandardAsyncJob, StandardAsyncRunInfo, StandardAsyncRunState]

  def mapInputPath(path: Path, tesJobPaths: TesJobPaths, commandDirectory: Path): String =
    path match {
      case drsPath: DrsPath =>
        val filepath = DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync()
        tesJobPaths.containerExec(commandDirectory, filepath)
      case httpPath: HttpPath =>
        // Strip the query params and anything after a # from http paths when turning them into local paths
        tesJobPaths.callInputsDockerRoot
          .resolve(httpPath.pathWithoutSchemeOrQueryOrFragment.stripPrefix("/"))
          .pathAsString
      case path: Path if path.startsWith(tesJobPaths.workflowPaths.DockerRoot) =>
        path.pathAsString
      case path: Path if path.equals(tesJobPaths.callExecutionRoot) =>
        commandDirectory.pathAsString
      case path: Path if path.startsWith(tesJobPaths.callExecutionRoot) =>
        tesJobPaths.containerExec(commandDirectory, path.name)
      case path: Path if path.startsWith(tesJobPaths.callRoot) =>
        tesJobPaths.callDockerRoot.resolve(path.name).pathAsString
      case _ =>
        tesJobPaths.callInputsDockerRoot.resolve(path.pathWithoutScheme.stripPrefix("/")).pathAsString
    }

  def getTaskEndTime(
    taskLogs: Future[Option[TaskLog]]
  )(implicit ec: ExecutionContext): Future[Option[String]] =
    taskLogs map { optTaskLog: Option[TaskLog] =>
      for {
        taskLog: TaskLog <- optTaskLog
        endTime: String <- taskLog.end_time
      } yield endTime
    }

  def getErrorSeq(taskLogs: Future[Option[TaskLog]])(implicit ec: ExecutionContext): Future[Option[Seq[String]]] =
    taskLogs.map(e => e.map(_.system_logs.getOrElse(Seq.empty[String])))

  def pollTesStatus(
    handle: StandardAsyncPendingExecutionHandle,
    fetchCostData: Boolean,
    fetchFullTaskViewFn: (StandardAsyncPendingExecutionHandle) => Future[Task],
    fetchMinimalTaskViewFn: (StandardAsyncPendingExecutionHandle) => Future[MinimalTaskView],
    getTesStatusFn: (Option[String], Option[TesVmCostData], String) => TesRunStatus,
    tellMetadataFn: (Map[String, Any]) => Unit,
    getErrorLogsFn: (StandardAsyncPendingExecutionHandle) => Future[Seq[String]]
  )(implicit ec: ExecutionContext): Future[TesRunStatus] =
    for {
      status <- queryStatusAndMaybeCostData(handle,
                                            fetchCostData,
                                            fetchFullTaskViewFn,
                                            fetchMinimalTaskViewFn,
                                            getTesStatusFn,
                                            tellMetadataFn
      )
      errorLog <- status match {
        case Error(_, _) | Failed(_, _) => getErrorLogsFn(handle)
        case _ => Future.successful(Seq.empty[String])
      }
      statusWithLog = status match {
        case Error(_, _) => Error(errorLog, handle.previousState.flatMap(_.costData))
        case Failed(_, _) => Failed(errorLog, handle.previousState.flatMap(_.costData))
        case _ => status
      }
    } yield statusWithLog

  def queryStatusAndMaybeCostData(
    handle: StandardAsyncPendingExecutionHandle,
    fetchCostData: Boolean,
    fetchFullTaskViewFn: (StandardAsyncPendingExecutionHandle) => Future[Task],
    fetchMinimalTaskViewFn: (StandardAsyncPendingExecutionHandle) => Future[MinimalTaskView],
    getTesStatusFn: (Option[String], Option[TesVmCostData], String) => TesRunStatus,
    tellMetadataFn: (Map[String, Any]) => Unit
  )(implicit ec: ExecutionContext): Future[TesRunStatus] =
    if (fetchCostData) {
      val task = fetchFullTaskViewFn(handle)
      task map { t =>
        val tesVmCostData = for {
          responseLogs <- t.logs
          startTime <- responseLogs.headOption.map(_.start_time)
          vmCost <- responseLogs.headOption.map(_.metadata.flatMap(_.get("vm_price_per_hour_usd")))
          // NB: End time is omitted here so we don't keep polling for it while the task runs. It will be acquired with a separate request when the task completes.
          tesVmCostData = TesVmCostData(startTime, None, vmCost)
        } yield tesVmCostData

        tesVmCostData match {
          case Some(v) =>
            val state = t.state
            v.startTime.foreach(s => tellMetadataFn(Map(CallMetadataKeys.VmStartTime -> s)))
            v.vmCost.foreach(v => tellMetadataFn(Map(CallMetadataKeys.VmCostPerHour -> v)))
            getTesStatusFn(state, tesVmCostData, handle.pendingJob.jobId)
          case None =>
            getTesStatusFn(t.state, tesVmCostData, handle.pendingJob.jobId)
        }
      }
    } else {
      val minimalTaskView = fetchMinimalTaskViewFn(handle)
      val previousCostData = handle.previousState.flatMap(c => c.costData)
      minimalTaskView map { t =>
        val state = t.state
        getTesStatusFn(Option(state), previousCostData, handle.pendingJob.jobId)
      }
    }

  def onTaskComplete(runStatus: TesRunStatus,
                     handle: StandardAsyncPendingExecutionHandle,
                     getTaskLogsFn: StandardAsyncPendingExecutionHandle => Future[Option[TaskLog]],
                     tellMetadataFn: Map[String, Any] => Unit,
                     logger: LoggingAdapter
  )(implicit ec: ExecutionContext): Unit = {
    val logs = getTaskLogsFn(handle)
    val taskEndTime = getTaskEndTime(logs)

    val errors = for {
      errors <- runStatus match {
        case Error(_, _) | Failed(_, _) => getErrorSeq(logs)
        case _ => Future.successful(Option(Seq.empty[String]))
      }
    } yield errors

    errors.onComplete {
      case Success(r) =>
        if (r.nonEmpty) {
          r.map(r => tellMetadataFn(Map(CallMetadataKeys.Failures -> r)))
        }
      case Failure(e) => logger.error(e.getMessage)
    }

    taskEndTime.onComplete {
      case Success(result) =>
        result.foreach(r => tellMetadataFn(Map(CallMetadataKeys.VmEndTime -> r)))
      case Failure(e) => logger.error(e.getMessage)
    }
  }
}

class TesAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
    extends BackendJobLifecycleActor
    with StandardAsyncExecutionActor
    with TesJobCachingActorHelper
    with CromwellInstrumentation {
  implicit val actorSystem = context.system
  implicit val materializer = ActorMaterializer()

  override type StandardAsyncRunInfo = Any

  override type StandardAsyncRunState = TesRunStatus

  def statusEquivalentTo(thiz: StandardAsyncRunState)(that: StandardAsyncRunState): Boolean = thiz == that
  override lazy val pollBackOff: SimpleExponentialBackoff = tesConfiguration.pollBackoff
  override lazy val executeOrRecoverBackOff: SimpleExponentialBackoff = tesConfiguration.executeOrRecoverBackoff

  private lazy val realDockerImageUsed: String =
    jobDescriptor.maybeCallCachingEligible.dockerHash.getOrElse(runtimeAttributes.dockerImage)
  override lazy val dockerImageUsed: Option[String] = Option(realDockerImageUsed)

  private val tesEndpoint = workflowDescriptor.workflowOptions.getOrElse("endpoint", tesConfiguration.endpointURL)

  override lazy val jobTag: String = jobDescriptor.key.tag

  private val outputMode = validate {
    OutputMode.withName(
      configurationDescriptor.backendConfig
        .getAs[String]("output-mode")
        .getOrElse("granular")
        .toUpperCase
    )
  }

  override def mapCommandLineWomFile(womFile: WomFile): WomFile =
    womFile.mapFile(value =>
      (getPath(value), asAdHocFile(womFile)) match {
        case (Success(path: Path), Some(adHocFile)) =>
          // Ad hoc files will be placed directly at the root ("/cromwell_root/ad_hoc_file.txt") unlike other input files
          // for which the full path is being propagated ("/cromwell_root/path/to/input_file.txt")
          tesJobPaths.containerExec(commandDirectory, adHocFile.alternativeName.getOrElse(path.name))
        case _ => mapCommandLineJobInputWomFile(womFile).value
      }
    )

  override def mapCommandLineJobInputWomFile(womFile: WomFile): WomFile =
    womFile.mapFile(value =>
      getPath(value) match {
        case Success(path: Path) =>
          TesAsyncBackendJobExecutionActor.mapInputPath(path, tesJobPaths, commandDirectory)
        case _ =>
          value
      }
    )

  override lazy val commandDirectory: Path =
    runtimeAttributes.dockerWorkingDir match {
      case Some(path) => DefaultPathBuilder.get(path)
      case None => tesJobPaths.callExecutionDockerRoot
    }

  def createTaskMessage(): ErrorOr[Task] = {
    val tesTask = (commandScriptContents, outputMode).mapN { case (contents, mode) =>
      TesTask(
        jobDescriptor,
        configurationDescriptor,
        jobLogger,
        tesJobPaths,
        runtimeAttributes,
        commandDirectory,
        contents,
        instantiatedCommand,
        realDockerImageUsed,
        mapCommandLineWomFile,
        jobShell,
        mode
      )
    }

    tesTask.map(TesTask.makeTask)
  }

  def writeScriptFile(): Future[Unit] =
    commandScriptContents.fold(
      errors => Future.failed(new RuntimeException(errors.toList.mkString(", "))),
      asyncIo.writeAsync(jobPaths.script, _, Seq.empty)
    )

  override def executeAsync(): Future[ExecutionHandle] = {
    // create call exec dir
    tesJobPaths.callExecutionRoot.createPermissionedDirectories()
    val taskMessageFuture = createTaskMessage().fold(
      errors => Future.failed(new RuntimeException(errors.toList.mkString(", "))),
      Future.successful
    )
    for {
      _ <- writeScriptFile()
      taskMessage <- taskMessageFuture
      entity <- Marshal(taskMessage).to[RequestEntity]
      ctr <- makeRequest[CreateTaskResponse](HttpRequest(method = HttpMethods.POST, uri = tesEndpoint, entity = entity))
    } yield PendingExecutionHandle(jobDescriptor, StandardAsyncJob(ctr.id), None, previousState = None)
  }

  override def reconnectAsync(jobId: StandardAsyncJob) = {
    val handle = PendingExecutionHandle[StandardAsyncJob, StandardAsyncRunInfo, StandardAsyncRunState](jobDescriptor,
                                                                                                       jobId,
                                                                                                       None,
                                                                                                       previousState =
                                                                                                         None
    )
    Future.successful(handle)
  }

  override def recoverAsync(jobId: StandardAsyncJob) = reconnectAsync(jobId)

  override def reconnectToAbortAsync(jobId: StandardAsyncJob) = {
    tryAbort(jobId)
    reconnectAsync(jobId)
  }

  override def tryAbort(job: StandardAsyncJob): Unit = {

    val returnCodeTmp = jobPaths.returnCode.plusExt("kill")
    returnCodeTmp.write(s"$SIGTERM\n")
    try {
      returnCodeTmp.copyTo(jobPaths.returnCode) // Will throw if file exists
      returnCodeTmp.delete(true)
    } catch {
      case _: FileAlreadyExistsException =>
        // If the process has already completed, there will be an existing rc file.
        returnCodeTmp.delete(true)
    }

    makeRequest[CancelTaskResponse](
      HttpRequest(method = HttpMethods.POST, uri = s"$tesEndpoint/${job.jobId}:cancel")
    ) onComplete {
      case Success(_) => jobLogger.info("{} Aborted {}", tag: Any, job.jobId)
      case Failure(ex) => jobLogger.warn("{} Failed to abort {}", tag, job.jobId, ex)
    }

    ()
  }

  override def requestsAbortAndDiesImmediately: Boolean = false

  override def onTaskComplete(runStatus: TesRunStatus, handle: StandardAsyncPendingExecutionHandle): Unit = {
    TesAsyncBackendJobExecutionActor.onTaskComplete(
      runStatus,
      handle,
      getTaskLogs,
      tellMetadata,
      log
    )
    ()
  }

  /*
   * We are polling for the status of a task to dynamically add its cost information to the metadata. pollStatusAsync
   * looks for a previous state on the handle; if it doesn't find anything (task just started), we poll for the status
   * and also fetch the cost data. If there is a previous status, we look and see if the cost data has been fetched.
   * If not, we poll for the status AND cost information in TES.
   * */
  override def pollStatusAsync(handle: StandardAsyncPendingExecutionHandle): Future[TesRunStatus] = {
    val fetchCostData = !handle.previousState.flatMap(_.costData).map(_.fullyPopulated).getOrElse(false)
    pollTesStatus(handle,
                  fetchCostData,
                  fetchFullTaskView,
                  fetchMinimalTesTask,
                  getTesStatus,
                  tellMetadata,
                  getErrorLogs
    )
  }

  private def getTesStatus(state: Option[String], withCostData: Option[TesVmCostData], jobId: String): TesRunStatus =
    state match {
      case s if s.contains("COMPLETE") =>
        jobLogger.info(s"Job ${jobId} is complete")
        Complete(withCostData)

      case s if s.contains("CANCELED") =>
        jobLogger.info(s"Job ${jobId} was canceled")
        Cancelled(withCostData)

      case s if s.contains("EXECUTOR_ERROR") =>
        jobLogger.info(s"TES reported a failure for Job ${jobId}: '$s'")
        Failed()

      case s if s.contains("SYSTEM_ERROR") =>
        jobLogger.info(s"TES reported an error for Job ${jobId}: '$s'")
        Error()

      case _ => Running(withCostData)
    }

  private def getErrorLogs(handle: StandardAsyncPendingExecutionHandle): Future[Seq[String]] = {
    val task = fetchFullTaskView(handle)
    task.map(t => t.logs.flatMap(_.lastOption).flatMap(_.system_logs).getOrElse(Seq.empty[String]))
  }

  private def getTaskLogs(handle: StandardAsyncPendingExecutionHandle): Future[Option[TaskLog]] = {
    val task = fetchFullTaskView(handle)
    val errorStates = List("EXECUTOR_ERROR", "SYSTEM_ERROR")

    val stateAndLogs = for {
      task <- task
      state = task.state
      taskLog = task.logs
    } yield (state, taskLog)

    stateAndLogs map { case (state, taskLog) =>
      if (errorStates.contains(state)) {
        Future.failed(new RuntimeException(s"Failed TES request: $state"))
      }
      taskLog.flatMap(_.headOption)
    }
  }

  private def fetchMinimalTesTask(handle: StandardAsyncPendingExecutionHandle): Future[MinimalTaskView] =
    makeRequest[MinimalTaskView](HttpRequest(uri = s"$tesEndpoint/${handle.pendingJob.jobId}?view=MINIMAL"))

  private def fetchFullTaskView(handle: StandardAsyncPendingExecutionHandle): Future[Task] =
    makeRequest[Task](HttpRequest(uri = s"$tesEndpoint/${handle.pendingJob.jobId}?view=FULL"))

  override def customPollStatusFailure: PartialFunction[(ExecutionHandle, Exception), ExecutionHandle] = {
    case (oldHandle: StandardAsyncPendingExecutionHandle @unchecked, e: Exception) =>
      jobLogger.error(s"$tag TES Job ${oldHandle.pendingJob.jobId} has not been found, failing call")
      FailedNonRetryableExecutionHandle(e, kvPairsToSave = None)
  }

  private def handleExecutionError(status: TesRunStatus, returnCode: Option[Int]): Future[ExecutionHandle] = {
    val msg = s"Task ${jobDescriptor.key.tag} failed with ${status.toString}"
    jobLogger.info(s"${msg}. Error messages: ${status.sysLogs}")
    val exception = new AggregatedMessageException(msg, status.sysLogs)
    Future.successful(FailedNonRetryableExecutionHandle(exception, returnCode, None))
  }

  override def handleExecutionFailure(status: StandardAsyncRunState, returnCode: Option[Int]) =
    status match {
      case Cancelled(_) => Future.successful(AbortedExecutionHandle)
      case Error(_, _) | Failed(_, _) => handleExecutionError(status, returnCode)
      case _ => super.handleExecutionFailure(status, returnCode)
    }

  override def isTerminal(runStatus: TesRunStatus): Boolean =
    runStatus.isTerminal

  override def isDone(runStatus: TesRunStatus): Boolean =
    runStatus match {
      case Complete(_) => true
      case _ => false
    }

  override def mapOutputWomFile(womFile: WomFile): WomFile =
    womFile mapFile { path =>
      val absPath = getPath(path) match {
        case Success(absoluteOutputPath) if absoluteOutputPath.isAbsolute => absoluteOutputPath
        case _ => tesJobPaths.callExecutionRoot.resolve(path)
      }

      if (!absPath.exists) {
        throw new FileNotFoundException(s"Could not process output, file not found: ${absPath.pathAsString}")
      } else absPath.pathAsString
    }

  // Headers that should be included with all requests to the TES server
  private def requestHeaders: List[HttpHeader] =
    tesConfiguration.token.flatMap { t =>
      HttpHeader.parse("Authorization", t) match {
        case Ok(header, _) => Some(header)
        case _ => None
      }
    }.toList

  private def makeRequest[A](request: HttpRequest)(implicit
    um: Unmarshaller[ResponseEntity, A]
  ): Future[A] =
    for {
      response <- withRetry(() => Http().singleRequest(request.withHeaders(requestHeaders)))
      data <-
        if (response.status.isFailure()) {
          response.entity.dataBytes.runFold(ByteString(""))(_ ++ _).map(_.utf8String) flatMap { errorBody =>
            Future.failed(
              new RuntimeException(s"Failed TES request: Code ${response.status.intValue()}, Body = $errorBody")
            )
          }
        } else {
          Unmarshal(response.entity).to[A]
        }
    } yield data

  override def platform: Option[Platform] = tesConfiguration.platform
}
