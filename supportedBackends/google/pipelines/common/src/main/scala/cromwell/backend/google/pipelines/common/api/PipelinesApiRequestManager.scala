package cromwell.backend.google.pipelines.common.api

import java.io.IOException
import java.util.UUID

import akka.actor.{Actor, ActorLogging, ActorRef, Props, SupervisorStrategy, Terminated, Timers}
import akka.dispatch.ControlMessage
import cats.data.NonEmptyList
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.{HttpHeaders, HttpRequest}
import common.util.Backoff
import cromwell.backend.BackendSingletonActorAbortWorkflow
import cromwell.backend.google.pipelines.common.PapiInstrumentation
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager._
import cromwell.backend.google.pipelines.common.api.clients.PipelinesApiRunCreationClient.JobAbortedException
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.core.{CromwellFatalExceptionMarker, LoadConfig, Mailbox, WorkflowId}
import cromwell.services.instrumentation.CromwellInstrumentationScheduler
import cromwell.services.loadcontroller.LoadControllerService.{HighLoad, LoadMetric, NormalLoad}
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric._

import scala.collection.immutable.Queue
import scala.concurrent.duration._
import scala.util.control.NoStackTrace

/**
  * Holds a set of PAPI request until a PipelinesApiRequestActor pulls the work.
  */
class PipelinesApiRequestManager(val qps: Int Refined Positive, requestWorkers: Int Refined Positive, override val serviceRegistryActor: ActorRef)
                                   (implicit batchHandler: PipelinesApiRequestHandler) extends Actor
  with ActorLogging with PapiInstrumentation with CromwellInstrumentationScheduler with Timers {

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
  */
  private val maxBatchRequestSize: Long = 14L * 1024L * 1024L
  private val requestTooLargeException = new UserPAPIApiException(
    cause = new IllegalArgumentException(s"The task run request has exceeded the maximum PAPI request size ($maxBatchRequestSize bytes)."),
    helpfulHint = Option("If you have a task with a very large number of inputs and / or outputs in your workflow you should try to reduce it. " +
      "Depending on your case you could: 1) Zip your input files together and unzip them in the command. 2) Use a file of file names " +
      "and localize the files yourself.")
  )

  private[api] lazy val nbWorkers = requestWorkers.value

  //
  private lazy val workerBatchInterval = determineBatchInterval(qps) * nbWorkers.toLong


  // workQueue is protected for the unit tests, not intended to be generally overridden
  protected[api] var workQueue: Queue[PAPIApiRequest] = Queue.empty
  private var workInProgress: Map[ActorRef, PipelinesApiWorkBatch] = Map.empty
  // Queries that have been scheduled to be retried through the actor's timer. They will boomerang back to this actor after
  // the scheduled delay, unless the workflow is aborted in the meantime in which case they will be cancelled.
  private var queriesWaitingForRetry: Set[PAPIApiRequest] = Set.empty[PAPIApiRequest]

  private def papiRequestWorkerProps = PipelinesApiRequestWorker.props(self, workerBatchInterval, serviceRegistryActor).withMailbox(Mailbox.PriorityMailbox)

  // statusPollers is protected for the unit tests, not intended to be generally overridden
  protected[api] var statusPollers: Vector[ActorRef] = resetAllWorkers()

  override def preStart() = {
    log.info("Running with {} PAPI request workers", requestWorkers.value)
    startInstrumentationTimer()
    super.preStart()
  }

  def monitorQueueSize() = {
    val load = if (workQueue.size > LoadConfig.PAPIThreshold) HighLoad else NormalLoad
    serviceRegistryActor ! LoadMetric("PAPIQueryManager", load)
    updateQueueSize(workQueue.size)
  }

  val requestManagerReceive: Receive = {
    case BackendSingletonActorAbortWorkflow(id) => abort(id)
    case status: PAPIStatusPollRequest => workQueue :+= status
    case create: PAPIRunCreationRequest =>
      if (create.contentLength > maxBatchRequestSize) {
        create.requester ! PipelinesApiRunCreationQueryFailed(create, requestTooLargeException)
      } else workQueue :+= create
    case abort: PAPIAbortRequest => workQueue :+= abort
    case PipelinesWorkerRequestWork(maxBatchSize) => handleWorkerAskingForWork(sender, maxBatchSize)
    case failure: PAPIApiRequestFailed =>
      handleQueryFailure(failure)
    case Terminated(actorRef) => onFailure(actorRef, new RuntimeException("PipelinesApiRequestHandler actor termination caught by manager") with NoStackTrace)
    case other => log.error(s"Unexpected message from {} to ${this.getClass.getSimpleName}: {}", sender().path.name, other)
  }

  override def receive = instrumentationReceive(monitorQueueSize _).orElse(requestManagerReceive)

  private def abort(workflowId: WorkflowId) = {
    def aborted(query: PAPIRunCreationRequest) = query.requester ! PipelinesApiRunCreationQueryFailed(query, JobAbortedException)

    workQueue = workQueue.filterNot({
      case query: PAPIRunCreationRequest if query.workflowId == workflowId =>
        aborted(query)
        true
      case _ => false
    })

    queriesWaitingForRetry = queriesWaitingForRetry.filterNot({
      case query: PAPIRunCreationRequest if query.workflowId == workflowId =>
        timers.cancel(query)
        queriesWaitingForRetry = queriesWaitingForRetry - query
        aborted(query)
        true
      case _ => false
    })
  }

  private def handleQueryFailure(failure: PAPIApiRequestFailed) = {
    val userError: Boolean = failure.cause.isInstanceOf[UserPAPIApiException]

    def failQuery(): Unit = {
      // NB: we don't count user errors towards the PAPI failed queries count in our metrics:
      if (!userError) {
        failedQuery(failure)
        log.warning(s"PAPI request workers tried and failed ${failure.query.failedAttempts} times to make ${failure.query.getClass.getSimpleName} request to PAPI")
      }

      failure.query.requester ! failure
    }

    def retryQuery(): Unit = {
      retriedQuery(failure)
      val nextRequest = failure.query.withFailedAttempt
      val delay = nextRequest.backoff.backoffMillis.millis
      queriesWaitingForRetry = queriesWaitingForRetry + nextRequest
      timers.startSingleTimer(nextRequest, nextRequest, delay)
    }

    if (userError || failure.query.failedAttempts >= maxRetries ) {
      failQuery()
    } else {
      retryQuery()
    }
  }

  private def handleWorkerAskingForWork(papiRequestWorkerActor: ActorRef, maxBatchSize: Int) = {
    log.debug("Request for PAPI requests received from {} (max batch size is {}, current queue size is {})", papiRequestWorkerActor.path.name, maxBatchSize, workQueue.size)

    workInProgress -= papiRequestWorkerActor
    val beheaded = beheadWorkQueue(maxBatchSize)
    beheaded.workToDo match {
      case Some(work) =>
        log.debug("Sending work to PAPI request worker {}", papiRequestWorkerActor.path.name)
        val workBatch = PipelinesApiWorkBatch(work)
        papiRequestWorkerActor ! workBatch
        workInProgress += (papiRequestWorkerActor -> workBatch)
      case None =>
        log.debug("No work for PAPI request worker {}", papiRequestWorkerActor.path.name)
        papiRequestWorkerActor ! NoWorkToDo
    }

    workQueue = beheaded.newWorkQueue
  }

  // Intentionally not final, this runs afoul of SI-4440 (I believe)
  private case class BeheadedWorkQueue(workToDo: Option[NonEmptyList[PAPIApiRequest]], newWorkQueue: Queue[PAPIApiRequest])
  private def beheadWorkQueue(maxBatchSize: Int): BeheadedWorkQueue = {
    import common.collections.EnhancedCollections._
    val DeQueued(head, tail) = workQueue.takeWhileWeighted(maxBatchRequestSize, _.contentLength, Option(maxBatchSize), strict = true)

    head.toList match {
      case h :: t => BeheadedWorkQueue(Option(NonEmptyList(h, t)), tail)
      case Nil => BeheadedWorkQueue(None, Queue.empty)
    }
  }

  /*
     Triggered by the 'case Terminated' receive handle whenever an actor being watched (ie a worker) terminates.
   */
  private def onFailure(terminee: ActorRef, throwable: => Throwable) = {
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
          case statusQuery: PAPIStatusPollRequest =>
            statusPolls += 1
            self ! PipelinesApiStatusQueryFailed(statusQuery, new SystemPAPIApiException(throwable))
          case runCreationQuery: PAPIRunCreationRequest =>
            runCreations += 1
            self ! PipelinesApiRunCreationQueryFailed(runCreationQuery, new SystemPAPIApiException(throwable))
          case abortQuery: PAPIAbortRequest =>
            abortRequests += 1
            self ! PipelinesApiAbortQueryFailed(abortQuery, new SystemPAPIApiException(throwable))
        }

        log.warning(
          s"PAPI request worker ${terminee.path.name} terminated. " +
            s"$runCreations run creation requests, $statusPolls status poll requests, and $abortRequests abort requests will be reconsidered. " +
            "If any of those succeeded in the cloud before the batch request failed, they might be run twice. " +
            s"Exception details: $throwable"
        )
      case None =>
        log.error(throwable, "The PAPI request worker '{}' terminated ({}). The request manager did not know what work the actor was doing so cannot resubmit any requests or inform any requesters of failures. This should never happen - please report this as a programmer error.", terminee.path.name, throwable.getMessage)
    }

    resetWorker(terminee)
  }

  private def resetWorker(worker: ActorRef): Unit = {

    if (!statusPollers.contains(worker)) {
      log.warning("PAPI request worker {} is being reset but was never registered in the pool of workers. This should never happen - please report this as a programmer error.", worker.path.name)
    }

    val stillGoodWorkers = statusPollers.filterNot(_ == worker)
    val newWorker = makeAndWatchWorkerActor()

    statusPollers = stillGoodWorkers :+ newWorker

    log.info("PAPI request worker {} has been removed and replaced by {} in the pool of {} workers", worker.path.name, newWorker.path.name, statusPollers.size)
  }

  private[api] def resetAllWorkers(): Vector[ActorRef] = {
    val pollers = Vector.fill(nbWorkers) { makeAndWatchWorkerActor() }
    pollers
  }

  private[api] def makeAndWatchWorkerActor(): ActorRef = {
    val newWorker = makeWorkerActor()
    context.watch(newWorker)
    newWorker
  }

  // Separate method to allow overriding in tests:
  private[api] def makeWorkerActor(): ActorRef = {
    context.actorOf(papiRequestWorkerProps, s"PAPIQueryWorker-${UUID.randomUUID()}")
  }
}

object PipelinesApiRequestManager {
  case object QueueMonitoringTimerKey
  case object QueueMonitoringTimerAction extends ControlMessage
  def props(qps: Int Refined Positive, requestWorkers: Int Refined Positive, serviceRegistryActor: ActorRef)
              (implicit batchHandler: PipelinesApiRequestHandler): Props = Props(new PipelinesApiRequestManager(qps, requestWorkers, serviceRegistryActor)).withDispatcher(BackendDispatcher)

  /**
    * Given the Genomics API queries per 100 seconds and given MaxBatchSize will determine a batch interval which
    * is at 90% of the quota. The (still crude) delta is to provide some room at the edges for things like new
    * calls, etc.
    */
  def determineBatchInterval(qps: Int Refined Positive): FiniteDuration = {
    val maxInterval = PipelinesApiRequestWorker.MaxBatchSize.toDouble / qps.value.toDouble
    val interval = ((maxInterval / 0.9) * 1000).toInt
    interval.milliseconds
  }

  sealed trait PAPIApiRequest {
    def workflowId: WorkflowId
    val failedAttempts: Int
    def requester: ActorRef
    def withFailedAttempt: PAPIApiRequest
    def backoff: Backoff
    def httpRequest: HttpRequest
    def contentLength: Long = (for {
      r <- Option(httpRequest)
      content <- Option(r.getContent)
    } yield content.getLength).getOrElse(0L)
  }
  private object PAPIApiRequest {
    // This must be a def, we want a new one each time (they're mutable! Boo!)
    def backoff: Backoff = SimpleExponentialBackoff(1.second, 1000.seconds, 1.5d)
  }

  case class PAPIStatusPollRequest(workflowId: WorkflowId,
                                   requester: ActorRef,
                                   httpRequest: HttpRequest,
                                   jobId: StandardAsyncJob,
                                   failedAttempts: Int = 0,
                                   backoff: Backoff = PAPIApiRequest.backoff) extends PAPIApiRequest {
    override def withFailedAttempt = this.copy(failedAttempts = failedAttempts + 1, backoff = backoff.next)
  }

  case class PAPIRunCreationRequest(workflowId: WorkflowId,
                                    requester: ActorRef,
                                    httpRequest: HttpRequest,
                                    failedAttempts: Int = 0,
                                    backoff: Backoff = PAPIApiRequest.backoff) extends PAPIApiRequest {
    override def withFailedAttempt = this.copy(failedAttempts = failedAttempts + 1, backoff = backoff.next)
  }

  case class PAPIAbortRequest(workflowId: WorkflowId,
                              requester: ActorRef,
                              httpRequest: HttpRequest,
                              jobId: StandardAsyncJob,
                              failedAttempts: Int = 0,
                              backoff: Backoff = PAPIApiRequest.backoff) extends PAPIApiRequest {
    override def withFailedAttempt = this.copy(failedAttempts = failedAttempts + 1, backoff = backoff.next)
  }

  sealed trait PAPIApiRequestFailed {
    val query: PAPIApiRequest
    val cause: PAPIApiException
  }

  final case class PipelinesApiStatusQueryFailed(query: PAPIApiRequest, cause: PAPIApiException) extends PAPIApiRequestFailed
  final case class PipelinesApiRunCreationQueryFailed(query: PAPIApiRequest, cause: PAPIApiException) extends PAPIApiRequestFailed
  final case class PipelinesApiAbortQueryFailed(query: PAPIApiRequest, cause: PAPIApiException) extends PAPIApiRequestFailed

  private[api] final case class PipelinesApiWorkBatch(workBatch: NonEmptyList[PAPIApiRequest])
  private[api] case object NoWorkToDo

  private[api] final case class PipelinesWorkerRequestWork(maxBatchSize: Int) extends ControlMessage

  final case class GoogleJsonException(e: GoogleJsonError, responseHeaders: HttpHeaders) extends IOException with CromwellFatalExceptionMarker {
    override def getMessage: String = e.getMessage
  }

  sealed trait PAPIApiException extends RuntimeException with CromwellFatalExceptionMarker {
    def cause: Throwable
    override def getCause: Throwable = cause
  }

  class SystemPAPIApiException(val cause: Throwable) extends PAPIApiException {
    override def getMessage: String = s"Unable to complete PAPI request due to system or connection error (${cause.getMessage})"
  }

  final class UserPAPIApiException(val cause: Throwable, helpfulHint: Option[String]) extends PAPIApiException {
    override def getMessage: String = s"Unable to complete PAPI request due to a problem with the request (${cause.getMessage}). ${helpfulHint.getOrElse("")}"
  }
}
