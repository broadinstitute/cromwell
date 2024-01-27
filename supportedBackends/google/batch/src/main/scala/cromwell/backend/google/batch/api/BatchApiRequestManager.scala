package cromwell.backend.google.batch.api

import java.io.IOException
import java.util.UUID
import akka.actor.{Actor, ActorLogging, ActorRef, Props, SupervisorStrategy, Terminated, Timers}
import akka.dispatch.ControlMessage
import cats.data.NonEmptyList
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpHeaders
import common.util.Backoff
import cromwell.backend.BackendSingletonActorAbortWorkflow
import cromwell.backend.google.batch.actors.BatchApiRunCreationClient.JobAbortedException
import cromwell.backend.google.batch.monitoring.BatchInstrumentation
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.core.{CromwellFatalExceptionMarker, LoadConfig, Mailbox, WorkflowId}
import cromwell.services.instrumentation.CromwellInstrumentationScheduler
import cromwell.services.loadcontroller.LoadControllerService.{HighLoad, LoadLevel, LoadMetric, NormalLoad}
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric._

import scala.collection.immutable.Queue
import scala.concurrent.duration._
import scala.util.control.NoStackTrace

// TODO: Alex - porting this file is almost done
/**
 * Holds a set of Batch request until a BatchApiRequestActor pulls the work.
 */
class BatchApiRequestManager(val qps: Int Refined Positive,
                             requestWorkers: Int Refined Positive,
                             override val serviceRegistryActor: ActorRef
)(implicit batchHandler: GcpBatchApiRequestHandler)
    extends Actor
    with ActorLogging
    with BatchInstrumentation
    with CromwellInstrumentationScheduler
    with Timers {

  import BatchApiRequestManager._

  override val supervisorStrategy = SupervisorStrategy.stoppingStrategy

  private val maxRetries = 10
  /*
   * Context: the batch.execute() method throws an IOException("insufficient data written") in certain conditions.
   * Here is what we know about it and how this attempts to address the issue.
   *
   * It was determined empirically that errors start to be thrown when the batch request approaches 15MB.
   * Looking more closely at timing it appears that the exception takes almost exactly 60 seconds to be thrown
   * from the batch.execute method, which suggests that this might be time related rather than byte size related and that
   * the 15MB limit is just an artifact of how much data can be sent / received in 60 seconds by the client / server.
   *
   * In an attempt to provide a fix for this issue, the total size of the batch size is limited to 14MB, which is a rather
   * arbitrary value only supported by local testing.
   *
   * Result of further investigation on the cause:
   * IOException("insufficient data written") is being thrown because the http request attempts to close() its output stream.
   * The close() method throws the "insufficient data written" exception because it still had data to send.
   * The close() method was called as part of a finally, because an exception was thrown earlier when attempting to write to the
   * stream. This exception is however swallowed by the one thrown in the close().
   * This commit https://github.com/google/google-http-java-client/pull/333 fixes the swallowing issue so the original
   * exception is thrown instead: IOException(“Error writing request body to server”).
   * Tracing back why this exception is being thrown, it appears that at some point the socket gets closed externally
   * (maybe the google server closes it ?)
   * which results in a SocketException("broken pipe") being thrown and eventually bubbles up to the IOExceptions above.
   *
   * see sun.net.www.protocol.http.HttpURLConnection
   * and com.google.api.client.http.javanet.NetHttpRequest
   *
   * TODO: Alex - do we still care about this? the github ticket was closed in 2016 and we don't use the batch API anyway
   */
  private val maxBatchRequestSize: Long = 14L * 1024L * 1024L // TODO: Alex what's the value for this?
  private val requestTooLargeException = new UserBatchApiException(
    cause = new IllegalArgumentException(
      s"The task run request has exceeded the maximum Batch request size ($maxBatchRequestSize bytes)."
    ),
    helpfulHint = Option(
      "If you have a task with a very large number of inputs and / or outputs in your workflow you should try to reduce it. " +
        "Depending on your case you could: 1) Zip your input files together and unzip them in the command. 2) Use a file of file names " +
        "and localize the files yourself."
    )
  )

  private[api] lazy val nbWorkers = requestWorkers.value

  private lazy val workerBatchInterval = determineBatchInterval(qps) * nbWorkers.toLong

  // workQueue is protected for the unit tests, not intended to be generally overridden
  protected[api] var workQueue: Queue[BatchApiRequest] = Queue.empty
  private var workInProgress: Map[ActorRef, BatchApiWorkBatch] = Map.empty
  // Queries that have been scheduled to be retried through the actor's timer. They will boomerang back to this actor after
  // the scheduled delay, unless the workflow is aborted in the meantime in which case they will be cancelled.
  private var queriesWaitingForRetry: Set[BatchApiRequest] = Set.empty[BatchApiRequest]

  private def batchRequestWorkerProps =
    BatchApiRequestWorker.props(self, workerBatchInterval, serviceRegistryActor).withMailbox(Mailbox.PriorityMailbox)

  // statusPollers is protected for the unit tests, not intended to be generally overridden
  protected[api] var statusPollers: Vector[ActorRef] = Vector.empty
  self ! ResetAllRequestWorkers

  private var previousLoad: LoadLevel = NormalLoad

  override def preStart() = {
    log.info("Running with {} Batch request workers", requestWorkers.value)
    startInstrumentationTimer()
    super.preStart()
  }

  private def monitorQueueSize(): Unit = {
    val newLoad = if (workQueue.size > LoadConfig.BatchThreshold) HighLoad else NormalLoad

    if (previousLoad == NormalLoad && newLoad == HighLoad)
      log.warning(
        s"Batch Request Manager transitioned to HighLoad with queue size ${workQueue.size} exceeding limit of ${LoadConfig.BatchThreshold}"
      )
    else if (previousLoad == HighLoad && newLoad == NormalLoad)
      log.info("Batch Request Manager transitioned back to NormaLoad")

    previousLoad = newLoad

    serviceRegistryActor ! LoadMetric("BatchQueryManager", newLoad)
    updateQueueSize(workQueue.size)
  }

  private val requestManagerReceive: Receive = {
    case ResetAllRequestWorkers => resetAllWorkers()
    case BackendSingletonActorAbortWorkflow(id) => abort(id)
    case status: BatchStatusPollRequest => workQueue :+= status
    case create: BatchRunCreationRequest =>
      if (create.contentLength > maxBatchRequestSize) {
        create.requester ! BatchApiRunCreationQueryFailed(create, requestTooLargeException)
      } else workQueue :+= create
    case abort: BatchAbortRequest => workQueue :+= abort
    case BatchWorkerRequestWork(maxBatchSize) => handleWorkerAskingForWork(sender(), maxBatchSize)
    case failure: BatchApiRequestFailed =>
      handleQueryFailure(failure)
    case Terminated(actorRef) =>
      onFailure(actorRef,
                new RuntimeException("BatchApiRequestHandler actor termination caught by manager") with NoStackTrace
      )
    case other =>
      log.error(s"Unexpected message from {} to ${this.getClass.getSimpleName}: {}", sender().path.name, other)
  }

  override def receive = instrumentationReceive(monitorQueueSize _).orElse(requestManagerReceive)

  private def abort(workflowId: WorkflowId): Unit = {
    def aborted(query: BatchRunCreationRequest) =
      query.requester ! BatchApiRunCreationQueryFailed(query, JobAbortedException)

    workQueue = workQueue.filterNot {
      case query: BatchRunCreationRequest if query.workflowId == workflowId =>
        aborted(query)
        true
      case _ => false
    }

    queriesWaitingForRetry = queriesWaitingForRetry.filterNot {
      case query: BatchRunCreationRequest if query.workflowId == workflowId =>
        timers.cancel(query)
        queriesWaitingForRetry = queriesWaitingForRetry - query
        aborted(query)
        true
      case _ => false
    }
  }

  private def handleQueryFailure(failure: BatchApiRequestFailed): Unit = {
    val userError: Boolean = failure.cause.isInstanceOf[UserBatchApiException]

    def failQuery(): Unit = {
      // NB: we don't count user errors towards the Batch failed queries count in our metrics:
      if (!userError) {
        // TODO: Alex - enable me
//        failedQuery(failure)
        log.warning(
          s"Batch request workers tried and failed ${failure.query.failedAttempts} times to make ${failure.query.getClass.getSimpleName} request to Batch"
        )
      }

      failure.query.requester ! failure
    }

    def retryQuery(): Unit = {
      // TODO: Alex - enable me
//      retriedQuery(failure)
      val nextRequest = failure.query.withFailedAttempt
      val delay = nextRequest.backoff.backoffMillis.millis
      queriesWaitingForRetry = queriesWaitingForRetry + nextRequest
      timers.startSingleTimer(nextRequest, nextRequest, delay)
    }

    if (userError || failure.query.failedAttempts >= maxRetries) {
      failQuery()
    } else {
      retryQuery()
    }
  }

  private def handleWorkerAskingForWork(batchRequestWorkerActor: ActorRef, maxBatchSize: Int): Unit = {
    log.debug(
      "Request for Batch requests received from {} (max batch size is {}, current queue size is {})",
      batchRequestWorkerActor.path.name,
      maxBatchSize,
      workQueue.size
    )

    workInProgress -= batchRequestWorkerActor
    val beheaded = beheadWorkQueue(maxBatchSize)
    beheaded.workToDo match {
      case Some(work) =>
        log.debug("Sending work to Batch request worker {}", batchRequestWorkerActor.path.name)
        val workBatch = BatchApiWorkBatch(work)
        batchRequestWorkerActor ! workBatch
        workInProgress += (batchRequestWorkerActor -> workBatch)
      case None =>
        log.debug("No work for Batch request worker {}", batchRequestWorkerActor.path.name)
        batchRequestWorkerActor ! NoWorkToDo
    }

    workQueue = beheaded.newWorkQueue
  }

  // Intentionally not final, this runs afoul of SI-4440 (I believe)
  private case class BeheadedWorkQueue(workToDo: Option[NonEmptyList[BatchApiRequest]],
                                       newWorkQueue: Queue[BatchApiRequest]
  )
  private def beheadWorkQueue(maxBatchSize: Int): BeheadedWorkQueue = {
    import common.collections.EnhancedCollections._

    val DeQueued(head, tail) =
      workQueue.takeWhileWeighted(maxBatchRequestSize, _.contentLength, Option(maxBatchSize), strict = true)

    head.toList match {
      case h :: t => BeheadedWorkQueue(Option(NonEmptyList(h, t)), tail)
      case Nil => BeheadedWorkQueue(None, Queue.empty)
    }
  }

  /*
     Triggered by the 'case Terminated' receive handle whenever an actor being watched (ie a worker) terminates.
   */
  private def onFailure(terminee: ActorRef, throwable: => Throwable): Unit = {
    // We assume this is a polling actor. Might change in a future update:
    workInProgress.get(terminee) match {
      case Some(work) =>
        // Most likely due to an unexpected HTTP error.
        // Do some book-keeping, log the error, and then push the work back onto the queue and keep going
        workInProgress -= terminee
        var runCreations = 0
        var statusPolls = 0
        var abortRequests = 0

        work.workBatch.toList.foreach {
          case statusQuery: BatchStatusPollRequest =>
            statusPolls += 1
            self ! BatchApiStatusQueryFailed(statusQuery, new SystemBatchApiException(throwable))
          case runCreationQuery: BatchRunCreationRequest =>
            runCreations += 1
            self ! BatchApiRunCreationQueryFailed(runCreationQuery, new SystemBatchApiException(throwable))
          case abortQuery: BatchAbortRequest =>
            abortRequests += 1
            self ! BatchApiAbortQueryFailed(abortQuery, new SystemBatchApiException(throwable))
        }

        log.warning(
          s"Batch request worker ${terminee.path.name} terminated. " +
            s"$runCreations run creation requests, $statusPolls status poll requests, and $abortRequests abort requests will be reconsidered. " +
            "If any of those succeeded in the cloud before the batch request failed, they might be run twice. " +
            s"Exception details: $throwable"
        )
      case None =>
        log.error(
          throwable,
          "The Batch request worker '{}' terminated ({}). The request manager did not know what work the actor was doing so cannot resubmit any requests or inform any requesters of failures. This should never happen - please report this as a programmer error.",
          terminee.path.name,
          throwable.getMessage
        )
    }

    resetWorker(terminee)
  }

  private def resetWorker(worker: ActorRef): Unit = {
    if (!statusPollers.contains(worker)) {
      log.warning(
        "Batch request worker {} is being reset but was never registered in the pool of workers. This should never happen - please report this as a programmer error.",
        worker.path.name
      )
    }

    val stillGoodWorkers = statusPollers.filterNot(_ == worker)
    val newWorker = makeAndWatchWorkerActor()

    statusPollers = stillGoodWorkers :+ newWorker

    log.info("Batch request worker {} has been removed and replaced by {} in the pool of {} workers",
             worker.path.name,
             newWorker.path.name,
             statusPollers.size
    )
  }

  private[api] def resetAllWorkers(): Unit = {
    log.info("'resetAllWorkers()' called to fill vector with {} new workers", nbWorkers)
    val result = Vector.fill(nbWorkers)(makeAndWatchWorkerActor())
    statusPollers = result
  }

  private[api] def makeAndWatchWorkerActor(): ActorRef = {
    val newWorker = makeWorkerActor()
    context.watch(newWorker)
    newWorker
  }

  // Separate method to allow overriding in tests:
  private[api] def makeWorkerActor(): ActorRef = {
    val result = context.actorOf(batchRequestWorkerProps, s"BatchQueryWorker-${UUID.randomUUID()}")
    log.info(
      s"Request manager ${self.path.name} created new Batch request worker ${result.path.name} with batch interval of ${workerBatchInterval}"
    )
    result
  }
}

object BatchApiRequestManager {
  case object ResetAllRequestWorkers
  // TODO: Alex - should we remove these unused items? they are not unsed in PAPIv2 either
  case object QueueMonitoringTimerKey
  case object QueueMonitoringTimerAction extends ControlMessage
  def props(qps: Int Refined Positive, requestWorkers: Int Refined Positive, serviceRegistryActor: ActorRef)(implicit
    batchHandler: GcpBatchApiRequestHandler
  ): Props =
    Props(new BatchApiRequestManager(qps, requestWorkers, serviceRegistryActor)).withDispatcher(BackendDispatcher)

  // TODO: Alex, what about Batch API?
  /**
   * Given the Genomics API queries per 100 seconds and given MaxBatchSize will determine a batch interval which
   * is at 90% of the quota. The (still crude) delta is to provide some room at the edges for things like new
   * calls, etc.
   */
  def determineBatchInterval(qps: Int Refined Positive): FiniteDuration = {
    val maxInterval = BatchApiRequestWorker.MaxBatchSize.toDouble / qps.value.toDouble
    val interval = ((maxInterval / 0.9) * 1000).toInt
    interval.milliseconds
  }

  sealed trait BatchApiRequest {
    def workflowId: WorkflowId
    val failedAttempts: Int
    def requester: ActorRef
    def withFailedAttempt: BatchApiRequest
    def backoff: Backoff
    def contentLength: Long
  }
  private object BatchApiRequest {
    // This must be a def, we want a new one each time (they're mutable! Boo!)
    def backoff: Backoff = SimpleExponentialBackoff(1.second, 1000.seconds, 1.5d)
  }

  case class BatchStatusPollRequest(workflowId: WorkflowId,
                                    requester: ActorRef,
                                    httpRequest: com.google.cloud.batch.v1.GetJobRequest,
                                    jobId: StandardAsyncJob,
                                    failedAttempts: Int = 0,
                                    backoff: Backoff = BatchApiRequest.backoff
  ) extends BatchApiRequest {
    override def withFailedAttempt = this.copy(failedAttempts = failedAttempts + 1, backoff = backoff.next)

    def contentLength: Long = (for {
      r <- Option(httpRequest)
      size <- Option(r.getSerializedSize)
    } yield size).getOrElse(0L)
  }

  case class BatchRunCreationRequest(workflowId: WorkflowId,
                                     requester: ActorRef,
                                     httpRequest: com.google.cloud.batch.v1.CreateJobRequest,
                                     failedAttempts: Int = 0,
                                     backoff: Backoff = BatchApiRequest.backoff
  ) extends BatchApiRequest {

    override def withFailedAttempt = this.copy(failedAttempts = failedAttempts + 1, backoff = backoff.next)

    def contentLength: Long = (for {
      r <- Option(httpRequest)
      size <- Option(r.getSerializedSize)
    } yield size).getOrElse(0L)
  }

  case class BatchAbortRequest(workflowId: WorkflowId,
                               requester: ActorRef,
                               httpRequest: com.google.cloud.batch.v1.DeleteJobRequest,
                               jobId: StandardAsyncJob,
                               failedAttempts: Int = 0,
                               backoff: Backoff = BatchApiRequest.backoff
  ) extends BatchApiRequest {
    override def withFailedAttempt = this.copy(failedAttempts = failedAttempts + 1, backoff = backoff.next)

    def contentLength: Long = (for {
      r <- Option(httpRequest)
      size <- Option(r.getSerializedSize)
    } yield size).getOrElse(0L)
  }

  sealed trait BatchApiRequestFailed {
    val query: BatchApiRequest
    val cause: BatchApiException
  }

  final case class BatchApiStatusQueryFailed(query: BatchApiRequest, cause: BatchApiException)
      extends BatchApiRequestFailed
  final case class BatchApiRunCreationQueryFailed(query: BatchApiRequest, cause: BatchApiException)
      extends BatchApiRequestFailed
  final case class BatchApiAbortQueryFailed(query: BatchApiRequest, cause: BatchApiException)
      extends BatchApiRequestFailed

  final private[api] case class BatchApiWorkBatch(workBatch: NonEmptyList[BatchApiRequest])
  private[api] case object NoWorkToDo

  final private[api] case class BatchWorkerRequestWork(maxBatchSize: Int) extends ControlMessage

  // TODO: Alex - is this necessary?
  final case class GoogleJsonException(e: GoogleJsonError, responseHeaders: HttpHeaders)
      extends IOException
      with CromwellFatalExceptionMarker {
    override def getMessage: String = e.getMessage
  }

  sealed trait BatchApiException extends RuntimeException with CromwellFatalExceptionMarker {
    def cause: Throwable
    override def getCause: Throwable = cause
  }

  class SystemBatchApiException(val cause: Throwable) extends BatchApiException {
    override def getMessage: String =
      s"Unable to complete Batch request due to system or connection error (${cause.getMessage})"
  }

  final class UserBatchApiException(val cause: Throwable, helpfulHint: Option[String]) extends BatchApiException {
    override def getMessage: String =
      s"Unable to complete Batch request due to a problem with the request (${cause.getMessage}). ${helpfulHint.getOrElse("")}"
  }
}
