package cromwell.backend.standard

import java.nio.file.Path

import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.event.LoggingReceive
import better.files.File
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, BackendJobExecutionResponse}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.async.AsyncBackendJobExecutionActor.{ExecutionMode, JobId, Recover}
import cromwell.backend.async.{AbortedExecutionHandle, AsyncBackendJobExecutionActor, ExecutionHandle, FailedNonRetryableExecutionHandle, PendingExecutionHandle, SuccessfulExecutionHandle}
import cromwell.backend.validation.{ContinueOnReturnCode, ContinueOnReturnCodeFlag}
import cromwell.backend.wdl.Command
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptor, BackendJobLifecycleActor}
import cromwell.services.keyvalue.KeyValueServiceActor._
import cromwell.services.metadata.CallMetadataKeys
import wdl4s._
import wdl4s.expression.{PureStandardLibraryFunctions, WdlFunctions}
import wdl4s.values.WdlValue

import scala.concurrent.{ExecutionContext, ExecutionContextExecutor, Future, Promise}
import scala.util.{Failure, Success, Try}

/**
  * An extension of the generic AsyncBackendJobExecutionActor providing a standard abstract implementation of an
  * asynchronous polling backend.
  *
  * Backends supported by this trait will all have common behavior. If a backend implementor wishes to provide a custom
  * behavior, one should instead implement the various methods in `AsyncBackendJobExecutionActor`.
  *
  * NOTE: Unlike the parent trait `AsyncBackendJobExecutionActor`, this trait is subject to even more frequent updates
  * as the common behavior among the backends adjusts in unison.
  */
trait StandardAsyncExecutionActor extends AsyncBackendJobExecutionActor {
  this: Actor with ActorLogging with BackendJobLifecycleActor =>

  val SIGTERM = 143
  val SIGINT = 130

  /** The type of the run info when a job is started. */
  type StandardAsyncRunInfo

  /** The type of the run status returned during each poll. */
  type StandardAsyncRunStatus

  /** The pending execution handle for each poll. */
  type StandardAsyncPendingExecutionHandle =
    PendingExecutionHandle[StandardAsyncJob, StandardAsyncRunInfo, StandardAsyncRunStatus]

  /** Standard set of parameters passed to the backend. */
  val standardParams: StandardAsyncExecutionActorParams

  override lazy val jobDescriptor: BackendJobDescriptor = standardParams.jobDescriptor

  override lazy val configurationDescriptor: BackendConfigurationDescriptor = standardParams.configurationDescriptor

  override lazy val completionPromise: Promise[BackendJobExecutionResponse] = standardParams.completionPromise

  /** Backend initialization data created by the a factory initializer. */
  lazy val backendInitializationDataOption: Option[BackendInitializationData] =
    standardParams.backendInitializationDataOption

  /** Typed backend initialization. */
  def backendInitializationDataAs[A <: BackendInitializationData]: A =
    BackendInitializationData.as[A](backendInitializationDataOption)

  /** @see [[StandardJobExecutionActorParams.serviceRegistryActor]] */
  lazy val serviceRegistryActor: ActorRef = standardParams.serviceRegistryActor

  /** @see [[StandardJobExecutionActorParams.jobIdKey]] */
  def jobIdKey: String = standardParams.jobIdKey

  /** @see [[Command.instantiate]] */
  def commandLineFunctions: WdlFunctions[WdlValue] = PureStandardLibraryFunctions

  /** @see [[Command.instantiate]] */
  def commandLinePreProcessor: EvaluatedTaskInputs => Try[EvaluatedTaskInputs] = Success.apply

  /** @see [[Command.instantiate]] */
  def commandLineValueMapper: WdlValue => WdlValue = identity

  /** The instantiated command. */
  lazy val instantiatedCommand: String = Command.instantiate(
    jobDescriptor, commandLineFunctions, commandLinePreProcessor, commandLineValueMapper).get

  /** A tag that may be used for logging. */
  lazy val tag = s"${this.getClass.getSimpleName} [UUID(${workflowId.shortString}):${jobDescriptor.key.tag}]"

  /**
    * When returns true, the `remoteStdErrPath` will be read. If contents of that path are non-empty, the job will fail.
    *
    * @return True if a non-empty `remoteStdErrPath` should fail the job.
    */
  def failOnStdErr: Boolean = false

  /**
    * Returns the path to the standard error output of the job. Only needs to be implemented if `failOnStdErr` is
    * returning `true`.
    *
    * @return The path to the standard error output.
    */
  def remoteStdErrPath: Path = {
    throw new NotImplementedError(s"failOnStdErr returned true but remote path not implemented by $getClass")
  }

  /**
    * Returns the path to the return code output of the job. Must be implemented unless `returnCodeContents` is
    * overridden not to use this method.
    *
    * @return The path to the return code output.
    */
  def remoteReturnCodePath: Path = {
    throw new NotImplementedError(s"remoteReturnCodePath returned true but remote path not implemented by $getClass")
  }

  /**
    * Returns the contents of the return code file.
    *
    * @return The contents of the return code file.
    */
  def returnCodeContents: String = File(remoteReturnCodePath).contentAsString

  /**
    * Returns the behavior for continuing on the return code, obtained by converting `returnCodeContents` to an Int.
    *
    * @return the behavior for continuing on the return code.
    */
  def continueOnReturnCode: ContinueOnReturnCode = ContinueOnReturnCodeFlag(false)

  /**
    * Returns the metadata key values to store before executing a job.
    *
    * @return the metadata key values to store before executing a job.
    */
  def startMetadataKeyValues: Map[String, Any] = Map.empty

  /**
    * Execute the job specified in the params. Should return a `StandardAsyncPendingExecutionHandle`, or a
    * `FailedExecutionHandle`.
    *
    * @return the execution handle for the job.
    */
  def execute(): ExecutionHandle

  /**
    * Recovers the specified job id, or starts a new job. The default implementation simply calls execute().
    *
    * @param jobId The previously recorded job id.
    * @return the execution handle for the job.
    */
  def recover(jobId: StandardAsyncJob): ExecutionHandle = execute()

  /**
    * Returns the run status for the job.
    *
    * @param handle The handle of the running job.
    * @return The status of the job.
    */
  def pollStatus(handle: StandardAsyncPendingExecutionHandle): StandardAsyncRunStatus = {
    throw new NotImplementedError(s"pollStatus nor pollStatusAsync not implemented by $getClass")
  }

  /**
    * Returns the async run status for the job.
    *
    * @param handle The handle of the running job.
    * @return The status of the job.
    */
  def pollStatusAsync(handle: StandardAsyncPendingExecutionHandle)
                     (implicit ec: ExecutionContext): Future[StandardAsyncRunStatus] = {
    Future.fromTry(Try(pollStatus(handle)))
  }

  /**
    * Adds custom behavior invoked when polling fails due to some exception. By default adds nothing.
    *
    * Examples may be when certain error codes are detected during polling, and a specific handle should be returned.
    *
    * For example, say a custom JobNotFoundException should be mapped to a FailedNonRetryableExecutionHandle.
    *
    * @return A partial function handler for exceptions after polling.
    */
  def customPollStatusFailure: PartialFunction[(ExecutionHandle, Exception), ExecutionHandle] = {
    PartialFunction.empty
  }

  /**
    * Maps the status to a status string that should be stored in the metadata.
    *
    * @param runStatus The run status.
    * @return Some() string that should be stored, or None if nothing should be stored in the metadata.
    */
  def statusString(runStatus: StandardAsyncRunStatus): Option[String] = None

  /**
    * Returns true when a job is complete, either successfully or unsuccessfully.
    *
    * @param runStatus The run status.
    * @return True if the job has completed.
    */
  def isTerminal(runStatus: StandardAsyncRunStatus): Boolean

  /**
    * Returns true if the status represents a success.
    *
    * @param runStatus The run status.
    * @return True if the job is a success.
    */
  def isSuccess(runStatus: StandardAsyncRunStatus): Boolean = true

  /**
    * Returns any custom metadata from the polled status.
    *
    * @param runStatus The run status.
    * @return The job metadata.
    */
  def getTerminalMetadata(runStatus: StandardAsyncRunStatus): Map[String, Any] = Map.empty

  /**
    * Attempts to abort a job when an abort signal is retrieved.
    *
    * If `abortAndDieImmediately` is true, then the actor will die immediately after this method returns.
    *
    * @param jobId The job to abort.
    */
  def tryAbort(jobId: StandardAsyncJob): Unit = {}

  /**
    * Returns true if when an abort signal is retrieved, the actor makes an attempt to abort and then immediately stops
    * itself _without_ polling for an aborted status.
    *
    * The default is false.
    *
    * @return true if actor should request an abort and then die immediately.
    */
  def requestsAbortAndDiesImmediately: Boolean = false

  /**
    * Return true if the return code is an abort code.
    *
    * By default, return codes `SIGINT` and `SIGTERM` return true.
    *
    * @param returnCode The return code.
    * @return True if the return code is for an abort.
    */
  def isAbort(returnCode: Int): Boolean = returnCode == SIGINT || returnCode == SIGTERM

  /**
    * Custom behavior to run after an abort signal is processed.
    *
    * By default handles the behavior of `abortAndDieImmediately`.
    */
  def postAbort(): Unit = {
    if (requestsAbortAndDiesImmediately) {
      context.parent ! AbortedResponse(jobDescriptor.key)
      context.stop(self)
    }
  }

  /**
    * Process a successful run, as defined by `isSuccess`.
    *
    * @param runStatus  The run status.
    * @param handle     The execution handle.
    * @param returnCode The return code.
    * @return The execution handle.
    */
  def handleExecutionSuccess(runStatus: StandardAsyncRunStatus, handle: StandardAsyncPendingExecutionHandle,
                             returnCode: Int): ExecutionHandle

  /**
    * Process an unsuccessful run, as defined by `isSuccess`.
    *
    * @param runStatus The run status.
    * @param handle    The execution handle.
    * @return The execution handle.
    */
  def handleExecutionFailure(runStatus: StandardAsyncRunStatus,
                             handle: StandardAsyncPendingExecutionHandle): ExecutionHandle = {
    FailedNonRetryableExecutionHandle(new Exception(s"Task failed for unknown reason: $runStatus"), None)
  }

  context.become(standardReceiveBehavior(None) orElse receive)

  def standardReceiveBehavior(jobIdOption: Option[StandardAsyncJob]): Receive = LoggingReceive {
    case AbortJobCommand =>
      jobIdOption foreach { jobId =>
        Try(tryAbort(jobId)) match {
          case Success(_) => jobLogger.info("{} Aborted {}", tag: Any, jobId)
          case Failure(ex) => jobLogger.warn("{} Failed to abort {}: {}", tag, jobId, ex.getMessage)
        }
      }
      postAbort()
    case KvPutSuccess(_) => // expected after the KvPut for the operation ID
  }

  override def retryable = false

  override def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    Future.fromTry(Try {
      val result = mode match {
        case Recover(jobId: StandardAsyncJob@unchecked) =>
          recover(jobId)
        case _ =>
          tellMetadata(startMetadataKeyValues)
          execute()
      }
      result match {
        case handle: PendingExecutionHandle[
          StandardAsyncJob@unchecked, StandardAsyncRunInfo@unchecked, StandardAsyncRunStatus@unchecked] =>
          tellKvJobId(handle.pendingJob)
          jobLogger.info(s"job id: ${handle.pendingJob.jobId}")
          tellMetadata(Map(CallMetadataKeys.JobId -> handle.pendingJob.jobId))
          context.become(standardReceiveBehavior(Option(handle.pendingJob)) orElse receive)
        case _ => /* ignore */
      }
      result
    } recoverWith {
      case exception: Exception =>
        jobLogger.error(s"Error attempting to $mode", exception)
        Failure(exception)
    })
  }

  override def poll(previous: ExecutionHandle)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    previous match {
      case handle: PendingExecutionHandle[
        StandardAsyncJob@unchecked, StandardAsyncRunInfo@unchecked, StandardAsyncRunStatus@unchecked] =>

        jobLogger.debug(s"$tag Polling Job ${handle.pendingJob}")
        pollStatusAsync(handle) map {
          backendRunStatus =>
            handlePollSuccess(handle, backendRunStatus)
        } recover {
          case throwable =>
            handlePollFailure(handle, throwable)
        }
      case successful: SuccessfulExecutionHandle => Future.successful(successful)
      case failed: FailedNonRetryableExecutionHandle => Future.successful(failed)
      case badHandle => Future.failed(new IllegalArgumentException(s"Unexpected execution handle: $badHandle"))
    }
  }

  /**
    * Process a poll success.
    *
    * @param oldHandle The previous execution status.
    * @param status The updated status.
    * @return The updated execution handle.
    */
  def handlePollSuccess(oldHandle: StandardAsyncPendingExecutionHandle,
                        status: StandardAsyncRunStatus): ExecutionHandle = {
    val previousStatus = oldHandle.previousStatus
    if (!(previousStatus contains status)) {
      /*
      If this is the first time checking the status, we log the transition as '-' to 'currentStatus'. Otherwise just use
      the state names.
       */
      val prevStateName = previousStatus.map(_.toString).getOrElse("-")
      jobLogger.info(s"$tag Status change from $prevStateName to $status")
      statusString(status) foreach { statusMetadata =>
        tellMetadata(Map(CallMetadataKeys.BackendStatus -> statusMetadata))
      }
    }

    status match {
      case _ if isTerminal(status) =>
        val metadata = getTerminalMetadata(status)
        tellMetadata(metadata)
        handleExecutionResult(status, oldHandle)
      case s => oldHandle.copy(previousStatus = Option(s)) // Copy the current handle with updated previous status.
    }
  }

  /**
    * Process a poll failure.
    *
    * @param oldHandle The previous execution handle.
    * @param throwable The cause of the polling failure.
    * @return The updated execution handle.
    */
  def handlePollFailure(oldHandle: StandardAsyncPendingExecutionHandle,
                        throwable: Throwable): ExecutionHandle = {
    throwable match {
      case exception: Exception =>
        val handler: PartialFunction[(ExecutionHandle, Exception), ExecutionHandle] =
          customPollStatusFailure orElse {
            case (handle: ExecutionHandle, exception: Exception) =>
              // Log exceptions and return the original handle to try again.
              jobLogger.warn(s"Caught exception, retrying", exception)
              handle
          }
        handler((oldHandle, exception))
      case error: Error => throw error // JVM-ending calamity.
      case _: Throwable =>
        // Someone has subclassed or instantiated Throwable directly. Kill the job. They should be using an Exception.
        FailedNonRetryableExecutionHandle(throwable)
    }
  }

  /**
    * Process an execution result.
    *
    * @param status The execution status.
    * @param oldHandle The previous execution handle.
    * @return The updated execution handle.
    */
  def handleExecutionResult(status: StandardAsyncRunStatus,
                            oldHandle: StandardAsyncPendingExecutionHandle): ExecutionHandle = {
    try {
      if (isSuccess(status)) {

        lazy val stderrLength: Long = File(remoteStdErrPath).size
        lazy val returnCode: Try[Int] = Try(returnCodeContents).map(_.trim.toInt)
        status match {
          case _ if failOnStdErr && stderrLength.intValue > 0 =>
            // returnCode will be None if it couldn't be downloaded/parsed, which will yield a null in the DB
            FailedNonRetryableExecutionHandle(new RuntimeException(
              s"execution failed: stderr has length $stderrLength"), returnCode.toOption)
          case _ if returnCode.isFailure =>
            val exception = returnCode.failed.get
            jobLogger.warn(s"could not download return code file, retrying", exception)
            // Return handle to try again.
            oldHandle
          case _ if returnCode.isFailure =>
            FailedNonRetryableExecutionHandle(new RuntimeException(
              s"execution failed: could not parse return code as integer", returnCode.failed.get))
          case _ if isAbort(returnCode.get) =>
            AbortedExecutionHandle
          case _ if !continueOnReturnCode.continueFor(returnCode.get) =>
            val message = s"Call ${jobDescriptor.key.tag}: return code was ${returnCode.get}"
            FailedNonRetryableExecutionHandle(new RuntimeException(message), returnCode.toOption)
          case _ =>
            handleExecutionSuccess(status, oldHandle, returnCode.get)
        }

      } else {
        handleExecutionFailure(status, oldHandle)
      }
    } catch {
      case e: Exception =>
        jobLogger.warn("Caught exception processing job result, retrying", e)
        // Return the original handle to try again.
        oldHandle
    }
  }

  /**
    * Send the job id of the running job to the key value store.
    *
    * @param runningJob The running job.
    */
  def tellKvJobId(runningJob: StandardAsyncJob): Unit = {
    val kvJobKey =
      KvJobKey(jobDescriptor.key.call.fullyQualifiedName, jobDescriptor.key.index, jobDescriptor.key.attempt)
    val scopedKey = ScopedKey(jobDescriptor.workflowDescriptor.id, kvJobKey, jobIdKey)
    val kvValue = Option(runningJob.jobId)
    val kvPair = KvPair(scopedKey, kvValue)
    val kvPut = KvPut(kvPair)
    serviceRegistryActor ! kvPut
  }

  /**
    * Sends metadata to the metadata store.
    *
    * @param metadataKeyValues Key/Values to store.
    */
  def tellMetadata(metadataKeyValues: Map[String, Any]): Unit = {
    import cromwell.services.metadata.MetadataService.implicits.MetadataAutoPutter
    serviceRegistryActor.putMetadata(jobDescriptor.workflowDescriptor.id, Option(jobDescriptor.key), metadataKeyValues)
  }

  override protected implicit def ec: ExecutionContextExecutor = context.dispatcher
}

/**
  * Implements the marker trait for a job id using a String.
  *
  * @param jobId The job id.
  */
case class StandardAsyncJob(jobId: String) extends JobId
