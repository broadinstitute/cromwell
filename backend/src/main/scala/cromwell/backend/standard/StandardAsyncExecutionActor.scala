package cromwell.backend.standard

import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.event.LoggingReceive
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, BackendJobExecutionResponse}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.async.AsyncBackendJobExecutionActor.{ExecutionMode, JobId, Recover}
import cromwell.backend.async.{AbortedExecutionHandle, AsyncBackendJobExecutionActor, ExecutionHandle, FailedNonRetryableExecutionHandle, PendingExecutionHandle, ReturnCodeIsNotAnInt, StderrNonEmpty, SuccessfulExecutionHandle, WrongReturnCode}
import cromwell.backend.validation._
import cromwell.backend.wdl.Command
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptor, BackendJobLifecycleActor}
import cromwell.services.io.AsyncIo
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
trait StandardAsyncExecutionActor extends AsyncBackendJobExecutionActor with StandardCachingActorHelper with AsyncIo {
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
  def execute(): Future[ExecutionHandle]

  /**
    * Recovers the specified job id, or starts a new job. The default implementation simply calls execute().
    *
    * @param jobId The previously recorded job id.
    * @return the execution handle for the job.
    */
  def recover(jobId: StandardAsyncJob): Future[ExecutionHandle] = execute()

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

  context.become(standardReceiveBehavior(None, abortReady = false) orElse receive)

  /**
    * @param abortReady will be true as soon as the execution handle is available and jobIdOption had a chance to be set to Some value.
    */
  def standardReceiveBehavior(jobIdOption: Option[StandardAsyncJob], abortReady: Boolean): Receive = LoggingReceive {
    case AbortJobCommand =>
      (jobIdOption, abortReady) match {
          // If we got a jobId, then we can abort so no need to check abortReady
        case (Some(jobId), _) =>
          Try(tryAbort(jobId)) match {
            case Success(_) => jobLogger.info("{} Aborted {}", tag: Any, jobId)
            case Failure(ex) => jobLogger.warn("{} Failed to abort {}: {}", tag, jobId, ex.getMessage)
          }
          postAbort()
          // If we don't have a job id but we're not ready to abort yet, re-submit the abort message
        case (None, false) => self ! AbortJobCommand
        case _ => postAbort()
      }
      jobIdOption foreach { jobId =>
        Try(tryAbort(jobId)) match {
          case Success(_) => jobLogger.info("{} Aborted {}", tag: Any, jobId)
          case Failure(ex) => jobLogger.warn("{} Failed to abort {}: {}", tag, jobId, ex.getMessage)
        }
      }
    case KvPutSuccess(_) => // expected after the KvPut for the operation ID
  }

  override def retryable = false

  override def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    val result = mode match {
      case Recover(jobId: StandardAsyncJob@unchecked) =>
        recover(jobId)
      case _ =>
        tellMetadata(startMetadataKeyValues)
        execute()
    }
    result onComplete {
      case Success(handle: PendingExecutionHandle[
        StandardAsyncJob@unchecked, StandardAsyncRunInfo@unchecked, StandardAsyncRunStatus@unchecked]) =>
        tellKvJobId(handle.pendingJob)
        jobLogger.info(s"job id: ${handle.pendingJob.jobId}")
        tellMetadata(Map(CallMetadataKeys.JobId -> handle.pendingJob.jobId))
        context.become(standardReceiveBehavior(Option(handle.pendingJob), abortReady = true) orElse receive)
      case Failure(exception) => 
        jobLogger.error(s"Error attempting to $mode: ${exception.getMessage}", exception)
        context.become(standardReceiveBehavior(None, abortReady = true) orElse receive)
      case _ =>
        context.become(standardReceiveBehavior(None, abortReady = true) orElse receive)
    }
    
    result
  }

  override def poll(previous: ExecutionHandle)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    previous match {
      case handle: PendingExecutionHandle[
        StandardAsyncJob@unchecked, StandardAsyncRunInfo@unchecked, StandardAsyncRunStatus@unchecked] =>

        jobLogger.debug(s"$tag Polling Job ${handle.pendingJob}")
        pollStatusAsync(handle) flatMap {
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
                        status: StandardAsyncRunStatus): Future[ExecutionHandle] = {
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
      case s => Future.successful(oldHandle.copy(previousStatus = Option(s))) // Copy the current handle with updated previous status.
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
                            oldHandle: StandardAsyncPendingExecutionHandle): Future[ExecutionHandle] = {

    if (isSuccess(status)) {
      val returnCode: Future[String] = read(jobPaths.returnCode)
      val stderrSize: Future[Long] = size(jobPaths.stderr)

      returnCode onFailure {
        case exception => jobLogger.warn(s"could not download return code file, retrying", exception)
      }

      stderrSize onFailure {
        case exception => jobLogger.warn(s"could not get stderr file size, retrying", exception)
      }

      val returnCodeAndStderrSize = for {
        rc <- returnCode
        size <- stderrSize
      } yield (size, rc, Try(rc.trim.toInt))

      returnCodeAndStderrSize map {
        case (_, rc, Failure(e)) =>
          FailedNonRetryableExecutionHandle(
            ReturnCodeIsNotAnInt(jobDescriptor.key.tag, rc, jobPaths.stderr)
          )
        case (length, _, rcAsInt) if failOnStdErr && length.intValue > 0 =>
          FailedNonRetryableExecutionHandle(
            StderrNonEmpty(jobDescriptor.key.tag, length, jobPaths.stderr), rcAsInt.toOption
          )
        case (_, _, Success(rcAsInt)) if isAbort(rcAsInt) =>
          AbortedExecutionHandle
        case (_, _, Success(rcAsInt)) if !continueOnReturnCode.continueFor(rcAsInt) =>
          FailedNonRetryableExecutionHandle(
            WrongReturnCode(jobDescriptor.key.tag, rcAsInt, jobPaths.stderr), Option(rcAsInt)
          )
        case (_, _, Success(rcAsInt)) =>
          handleExecutionSuccess(status, oldHandle, rcAsInt)
      } recover {
        case e => oldHandle
      }
      
    } else Future.successful(handleExecutionFailure(status, oldHandle))
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
