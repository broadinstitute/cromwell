package cromwell.backend.google.pipelines.common.api

import java.io.IOException
import java.util.UUID

import akka.actor.{Actor, ActorLogging, ActorRef, Props, Terminated, Timers}
import akka.dispatch.ControlMessage
import cats.data.NonEmptyList
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.{HttpHeaders, HttpRequest}
import common.util.Backoff
import cromwell.backend.BackendSingletonActorAbortWorkflow
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager._
import cromwell.backend.google.pipelines.common.api.clients.PipelinesApiRunCreationClient.JobAbortedException
import cromwell.backend.google.pipelines.common.PapiInstrumentation
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.core.{CromwellFatalExceptionMarker, LoadConfig, Mailbox, WorkflowId}
import cromwell.services.instrumentation.{CromwellInstrumentation, CromwellInstrumentationScheduler}
import cromwell.services.loadcontroller.LoadControllerService.{HighLoad, LoadMetric, NormalLoad}
import cromwell.util.StopAndLogSupervisor
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric._

import scala.collection.immutable.Queue
import scala.concurrent.duration._

/**
  * Holds a set of JES API requests until a JesQueryActor pulls the work.
  */
class PipelinesApiRequestManager(val qps: Int Refined Positive, requestWorkers: Int Refined Positive, override val serviceRegistryActor: ActorRef)
                                   (implicit batchHandler: PipelinesApiRequestHandler) extends Actor
  with ActorLogging with StopAndLogSupervisor with PapiInstrumentation with CromwellInstrumentationScheduler with Timers {

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
  private val requestTooLargeException = new PAPIApiException(new IllegalArgumentException(
    "The task run request has exceeded the maximum PAPI request size." +
      "If you have a task with a very large number of inputs and / or outputs in your workflow you should try to reduce it. " +
      "Depending on your case you could: 1) Zip your input files together and unzip them in the command. 2) Use a file of file names " +
      "and localize the files yourself."
  ))

  private[api] lazy val nbWorkers = requestWorkers.value

  // 
  private lazy val workerBatchInterval = determineBatchInterval(qps) * nbWorkers.toLong

  scheduleInstrumentation { updateQueueSize(workQueue.size) }

  // workQueue is protected for the unit tests, not intended to be generally overridden
  protected[api] var workQueue: Queue[PAPIApiRequest] = Queue.empty
  private var workInProgress: Map[ActorRef, PipelinesApiWorkBatch] = Map.empty
  // Queries that have been scheduled to be retried through the actor's timer. They will boomerang back to this actor after
  // the scheduled delay, unless the workflow is aborted in the meantime in which case they will be cancelled.
  private var queriesWaitingForRetry: Set[PAPIApiRequest] = Set.empty[PAPIApiRequest]

  private def statusPollerProps = PipelinesApiRequestWorker.props(self, workerBatchInterval, serviceRegistryActor).withMailbox(Mailbox.PriorityMailbox)

  // statusPollers is protected for the unit tests, not intended to be generally overridden
  protected[api] var statusPollers: Vector[ActorRef] = resetAllWorkers()

  override def preStart() = {
    log.info("{} Running with {} workers", self.path.name, requestWorkers.value)
    timers.startSingleTimer(QueueMonitoringTimerKey, QueueMonitoringTimerAction, CromwellInstrumentation.InstrumentationRate)
    super.preStart()
  }

  def monitorQueueSize() = {
    val load = if (workQueue.size > LoadConfig.PAPIThreshold) HighLoad else NormalLoad
    serviceRegistryActor ! LoadMetric("PAPIQueryManager", load)
    timers.startSingleTimer(QueueMonitoringTimerKey, QueueMonitoringTimerAction, CromwellInstrumentation.InstrumentationRate)
  }

  override def receive = {
    case QueueMonitoringTimerAction => monitorQueueSize()
    case BackendSingletonActorAbortWorkflow(id) => abort(id)
    case status: PAPIStatusPollRequest => workQueue :+= status
    case create: PAPIRunCreationRequest =>
      if (create.contentLength > maxBatchRequestSize) {
        create.requester ! PipelinesApiRunCreationQueryFailed(create, requestTooLargeException)
      } else workQueue :+= create
    case abort: PAPIAbortRequest => workQueue :+= abort
    case PipelinesWorkerRequestWork(maxBatchSize) =>
      log.debug("Request for JES Polling Work received (max batch: {}, current queue size is {})", maxBatchSize, workQueue.size)
      handleJesPollingRequest(sender, maxBatchSize)
    case failure: PAPIApiRequestFailed => handleQueryFailure(failure)
    case Terminated(actorRef) => onFailure(actorRef, new RuntimeException("Polling stopped itself unexpectedly"))
    case other => log.error(s"Unexpected message to JesPollingManager: $other")
  }

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

  private def handleQueryFailure(failure: PAPIApiRequestFailed) = if (failure.query.failedAttempts < maxRetries) {
    val nextRequest = failure.query.withFailedAttempt
    val delay = nextRequest.backoff.backoffMillis.millis
    queriesWaitingForRetry = queriesWaitingForRetry + nextRequest
    timers.startSingleTimer(nextRequest, nextRequest, delay)
    failedQuery(failure)
  } else {
    retriedQuery(failure)
    failure.query.requester ! failure
  }

  private def handleJesPollingRequest(workPullingJesPollingActor: ActorRef, maxBatchSize: Int) = {
    workInProgress -= workPullingJesPollingActor
    val beheaded = beheadWorkQueue(maxBatchSize)
    beheaded.workToDo match {
      case Some(work) =>
        log.debug("Sending work to JES API query manager.")
        val workBatch = PipelinesApiWorkBatch(work)
        workPullingJesPollingActor ! workBatch
        workInProgress += (workPullingJesPollingActor -> workBatch)
      case None =>
        log.debug("No work for JES poller")
        workPullingJesPollingActor ! NoWorkToDo
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

  override protected def onFailure(terminee: ActorRef, throwable: => Throwable) = {
    // We assume this is a polling actor. Might change in a future update:
    workInProgress.get(terminee) match {
      case Some(work) =>
        // Most likely due to an unexpected HTTP error, push the work back on the queue and keep going
        log.error(throwable, s"The JES API worker actor $terminee unexpectedly terminated while conducting ${work.workBatch.tail.size + 1} polls. Making a new one...")
        workInProgress -= terminee
        work.workBatch.toList.foreach {
          case statusQuery: PAPIStatusPollRequest =>
            self ! PipelinesApiStatusQueryFailed(statusQuery, new PAPIApiException(throwable))
          case runCreationQuery: PAPIRunCreationRequest =>
            self ! PipelinesApiRunCreationQueryFailed(runCreationQuery, new PAPIApiException(throwable))
          case abortQuery: PAPIAbortRequest =>
            self ! PipelinesApiAbortQueryFailed(abortQuery, new PAPIApiException(throwable))
        }
      case None =>
        // It managed to die while doing absolutely nothing...!?
        // Maybe it deserves an entry in https://en.wikipedia.org/wiki/List_of_unusual_deaths
        // Oh well, in the mean time don't do anything, just start a new one
        log.error(throwable, s"The JES API worker actor managed to unexpectedly terminate whilst doing absolutely nothing (${throwable.getMessage}). This is probably a programming error. Making a new one...")
    }

    resetWorker(terminee)
  }

  private def resetWorker(worker: ActorRef) = {
    statusPollers = statusPollers.filterNot(_ == worker) :+ makeWorkerActor()
  }

  private[api] def resetAllWorkers(): Vector[ActorRef] = {
    val pollers = Vector.fill(nbWorkers) { makeWorkerActor() }
    pollers
  }

  // Separate method to allow overriding in tests:
  private[api] def makeWorkerActor(): ActorRef = {
    val newWorker = context.actorOf(statusPollerProps, s"PAPIQueryWorker-${UUID.randomUUID()}")
    context.watch(newWorker)
    newWorker
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
    def contentLength: Long = Option(httpRequest.getContent).map(_.getLength).getOrElse(0L)
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

  class PAPIApiException(val e: Throwable) extends RuntimeException(e) with CromwellFatalExceptionMarker {
    override def getMessage: String = "Unable to complete JES Api Request"
  }
}
