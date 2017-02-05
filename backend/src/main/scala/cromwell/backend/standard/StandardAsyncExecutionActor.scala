package cromwell.backend.standard

import java.nio.file.Path

import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.event.LoggingReceive
import better.files.File
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, BackendJobExecutionResponse}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.async.AsyncBackendJobExecutionActor.{ExecutionMode, JobId, Recover}
import cromwell.backend.async.{AbortedExecutionHandle, AsyncBackendJobExecutionActor, ExecutionHandle, FailedNonRetryableExecutionHandle, PendingExecutionHandle, ReturnCodeIsNotAnInt, StderrNonEmpty, SuccessfulExecutionHandle, WrongReturnCode}
import cromwell.backend.validation._
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
trait StandardAsyncExecutionActor extends AsyncBackendJobExecutionActor with StandardCachingActorHelper {
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
  def standardParams: StandardAsyncExecutionActorParams

  override lazy val configurationDescriptor: BackendConfigurationDescriptor = standardParams.configurationDescriptor

  override lazy val completionPromise: Promise[BackendJobExecutionResponse] = standardParams.completionPromise

  /** Backend initialization data created by the a factory initializer. */
  override lazy val backendInitializationDataOption: Option[BackendInitializationData] =
    standardParams.backendInitializationDataOption

  /** @see [[StandardAsyncExecutionActorParams.serviceRegistryActor]] */
  override lazy val serviceRegistryActor: ActorRef = standardParams.serviceRegistryActor

  /** @see [[StandardAsyncExecutionActorParams.jobDescriptor]] */
  override lazy val jobDescriptor: BackendJobDescriptor = standardParams.jobDescriptor

  /** @see [[StandardAsyncExecutionActorParams.jobIdKey]] */
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
  lazy val failOnStdErr: Boolean = RuntimeAttributesValidation.extract(
    FailOnStderrValidation.instance, validatedRuntimeAttributes)

  /**
    * Returns the behavior for continuing on the return code, obtained by converting `returnCodeContents` to an Int.
    *
    * @return the behavior for continuing on the return code.
    */
  lazy val continueOnReturnCode: ContinueOnReturnCode = RuntimeAttributesValidation.extract(
    ContinueOnReturnCodeValidation.instance, validatedRuntimeAttributes)

  /**
    * Execute the job specified in the params. Should return a `StandardAsyncPendingExecutionHandle`, or a
    * `FailedExecutionHandle`.
    *
    * @return the execution handle for the job.
    */
  def execute(): ExecutionHandle = {
    throw new NotImplementedError(s"Neither execute() nor executeAsync() implemented by $getClass")
  }

  /**
    * Async execute the job specified in the params. Should return a `StandardAsyncPendingExecutionHandle`, or a
    * `FailedExecutionHandle`.
    *
    * @return the execution handle for the job.
    */
  def executeAsync()(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    Future.fromTry(Try(execute()))
  }

  /**
    * Recovers the specified job id, or starts a new job. The default implementation simply calls execute().
    *
    * @param jobId The previously recorded job id.
    * @return the execution handle for the job.
    */
  def recover(jobId: StandardAsyncJob): ExecutionHandle = execute()

  /**
    * Async recovers the specified job id, or starts a new job. The default implementation simply calls execute().
    *
    * @param jobId The previously recorded job id.
    * @return the execution handle for the job.
    */
  def recoverAsync(jobId: StandardAsyncJob)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    Future.fromTry(Try(recover(jobId)))
  }

  /**
    * Returns the run status for the job.
    *
    * @param handle The handle of the running job.
    * @return The status of the job.
    */
  def pollStatus(handle: StandardAsyncPendingExecutionHandle): StandardAsyncRunStatus = {
    throw new NotImplementedError(s"Neither pollStatus nor pollStatusAsync implemented by $getClass")
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
                             handle: StandardAsyncPendingExecutionHandle,
                             returnCode: Option[Int]): ExecutionHandle = {
    FailedNonRetryableExecutionHandle(new Exception(s"Task failed for unknown reason: $runStatus"), returnCode)
  }

  // See executeOrRecoverSuccess
  private var missedAbort = false
  private case object CheckMissedAbort

  context.become(standardReceiveBehavior(None) orElse receive)

  def standardReceiveBehavior(jobIdOption: Option[StandardAsyncJob]): Receive = LoggingReceive {
    case AbortJobCommand =>
      jobIdOption match {
        case Some(jobId) =>
          Try(tryAbort(jobId)) match {
            case Success(_) => jobLogger.info("{} Aborted {}", tag: Any, jobId)
            case Failure(ex) => jobLogger.warn("{} Failed to abort {}: {}", tag, jobId, ex.getMessage)
          }
        case None => missedAbort = true
      }
      postAbort()
    case CheckMissedAbort =>
      if (missedAbort)
        self ! AbortJobCommand
    case KvPutSuccess(_) => // expected after the KvPut for the operation ID
  }

  override def retryable = false

  override def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    val executeOrRecoverFuture = {
      mode match {
        case Recover(jobId: StandardAsyncJob@unchecked) => recoverAsync(jobId)
        case _ =>
          tellMetadata(startMetadataKeyValues)
          executeAsync()
      }
    }

    executeOrRecoverFuture map executeOrRecoverSuccess recoverWith {
      case throwable: Throwable => Future failed {
        jobLogger.error(s"Error attempting to $mode", throwable)
        throwable
      }
    }
  }

  private def executeOrRecoverSuccess(executionHandle: ExecutionHandle): ExecutionHandle = {
    executionHandle match {
      case handle: PendingExecutionHandle[
        StandardAsyncJob@unchecked, StandardAsyncRunInfo@unchecked, StandardAsyncRunStatus@unchecked] =>
        tellKvJobId(handle.pendingJob)
        jobLogger.info(s"job id: ${handle.pendingJob.jobId}")
        tellMetadata(Map(CallMetadataKeys.JobId -> handle.pendingJob.jobId))
        /*
        NOTE: Because of the async nature of the Scala Futures, there is a point in time where we have submitted this or
        the prior runnable to the thread pool this actor doesn't know the job id for aborting. These runnables are
        queued up and may still be run by the thread pool anytime in the future. Issue #1218 may address this
        inconsistency at a later time. For now, just go back and check if we missed the abort command.
        */
        context.become(standardReceiveBehavior(Option(handle.pendingJob)) orElse receive)
        self ! CheckMissedAbort
      case _ =>
    }
    executionHandle
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
      tellMetadata(Map(CallMetadataKeys.BackendStatus -> status))
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

      lazy val returnCodeAsString: Try[String] = Try(File(jobPaths.returnCode).contentAsString)
      lazy val returnCodeAsInt: Try[Int] = returnCodeAsString.map(_.trim.toInt)
      lazy val stderrAsOption: Option[Path] = Option(jobPaths.stderr)
      
      if (isSuccess(status)) {
        lazy val stderrLength: Try[Long] = Try(File(jobPaths.stderr).size)
        (stderrLength, returnCodeAsString, returnCodeAsInt) match {
            // Failed to get stderr size -> Retry
          case (Failure(exception), _, _) =>
            jobLogger.warn(s"could not get stderr file size, retrying", exception)
            oldHandle
            // Failed to get return code content -> Retry
          case (_, Failure(exception), _) =>
            jobLogger.warn(s"could not download return code file, retrying", exception)
            oldHandle
            // Failed to convert return code content to Int -> Fail
          case (_, _, Failure(_)) =>
            FailedNonRetryableExecutionHandle(ReturnCodeIsNotAnInt(jobDescriptor.key.tag, returnCodeAsString.get, stderrAsOption))
            // Stderr is not empty and failOnStdErr is true -> Fail
          case (Success(length), _, _) if failOnStdErr && length.intValue > 0 =>
            FailedNonRetryableExecutionHandle(StderrNonEmpty(jobDescriptor.key.tag, length, stderrAsOption), returnCodeAsInt.toOption)
            // Return code is abort code -> Abort
          case (_, _, Success(rc)) if isAbort(rc) =>
            AbortedExecutionHandle
            // Return code is not valid -> Fail
          case (_, _, Success(rc)) if !continueOnReturnCode.continueFor(rc) =>
            FailedNonRetryableExecutionHandle(WrongReturnCode(jobDescriptor.key.tag, returnCodeAsInt.get, stderrAsOption), returnCodeAsInt.toOption)
            // Otherwise -> Succeed
          case (_, _, Success(rc)) => 
            handleExecutionSuccess(status, oldHandle, rc)
        }
      } else {
        handleExecutionFailure(status, oldHandle, returnCodeAsInt.toOption)
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
