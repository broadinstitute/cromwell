package cromwell.backend.impl.jes.statuspolling

import java.io.IOException

import akka.actor.{Actor, ActorLogging, ActorRef, Props, Terminated, Timers}
import cats.data.NonEmptyList
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.{HttpHeaders, HttpRequest}
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model.RunPipelineRequest
import cromwell.backend.BackendSingletonActorAbortWorkflow
import cromwell.backend.impl.jes.Run
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager.{JesApiException, _}
import cromwell.backend.impl.jes.statuspolling.JesRunCreationClient.JobAbortedException
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core.retry.{Backoff, SimpleExponentialBackoff}
import cromwell.core.{CromwellFatalExceptionMarker, WorkflowId}
import cromwell.services.instrumentation.CromwellInstrumentationScheduler
import cromwell.util.StopAndLogSupervisor
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric._

import scala.collection.immutable.Queue
import scala.concurrent.duration._

/**
  * Holds a set of JES API requests until a JesQueryActor pulls the work.
  */
class JesApiQueryManager(val qps: Int Refined Positive, override val serviceRegistryActor: ActorRef) extends Actor 
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
  private val requestTooLargeException = new JesApiException(new IllegalArgumentException(
    "The task run request has exceeded the maximum PAPI request size." +
      "If you have a task with a very large number of inputs and / or outputs in your workflow you should try to reduce it. " +
      "Depending on your case you could: 1) Zip your input files together and unzip them in the command. 2) Use a file of file names " +
      "and localize the files yourself."
  ))

  scheduleInstrumentation { updateQueueSize(workQueue.size) }

  // workQueue is protected for the unit tests, not intended to be generally overridden
  protected[statuspolling] var workQueue: Queue[JesApiQuery] = Queue.empty
  private var workInProgress: Map[ActorRef, JesPollingWorkBatch] = Map.empty
  // Queries that have been scheduled to be retried through the actor's timer. They will boomerang back to this actor after
  // the scheduled delay, unless the workflow is aborted in the meantime in which case they will be cancelled.
  private var queriesWaitingForRetry: Set[JesApiQuery] = Set.empty[JesApiQuery]

  private def statusPollerProps = JesPollingActor.props(self, qps, serviceRegistryActor)

  // statusPoller is protected for the unit tests, not intended to be generally overridden
  protected[statuspolling] var statusPoller: ActorRef = _

  resetWorker()

  override def receive = {
    case BackendSingletonActorAbortWorkflow(id) => abort(id)
    case DoPoll(workflowId, run) => workQueue :+= makePollQuery(workflowId, sender, run)
    case DoCreateRun(workflowId, genomics, rpr) =>
      val creationQuery = makeCreateQuery(workflowId, sender, genomics, rpr)

      if (creationQuery.contentLength > maxBatchRequestSize) {
        creationQuery.requester ! JesApiRunCreationQueryFailed(creationQuery, requestTooLargeException)
      } else workQueue :+= creationQuery
    case q: JesApiQuery => workQueue :+= q
    case RequestJesPollingWork(maxBatchSize) =>
      log.debug("Request for JES Polling Work received (max batch: {}, current queue size is {})", maxBatchSize, workQueue.size)
      handleJesPollingRequest(sender, maxBatchSize)
    case failure: JesApiQueryFailed => handleQueryFailure(failure)
    case Terminated(actorRef) => onFailure(actorRef, new RuntimeException("Polling stopped itself unexpectedly"))
    case other => log.error(s"Unexpected message to JesPollingManager: $other")
  }
  
  private def abort(workflowId: WorkflowId) = {
    def aborted(query: JesApiQuery) = query match {
      case creation: JesRunCreationQuery => creation.requester ! JesApiRunCreationQueryFailed(creation, JobAbortedException)
      case poll: JesStatusPollQuery => poll.requester ! JesApiStatusQueryFailed(poll, JobAbortedException)
    }

    workQueue = workQueue.filterNot({
      case query if query.workflowId == workflowId =>
        aborted(query)
        true
      case _ => false
    })

    queriesWaitingForRetry = queriesWaitingForRetry.filterNot({
      case query if query.workflowId == workflowId =>
        timers.cancel(query)
        queriesWaitingForRetry = queriesWaitingForRetry - query
        aborted(query)
        true
      case _ => false
    })
  }

  private [statuspolling] def makeCreateQuery(workflowId: WorkflowId, replyTo: ActorRef, genomics: Genomics, rpr: RunPipelineRequest) = {
    JesRunCreationQuery(workflowId, replyTo, genomics, rpr)
  }

  private [statuspolling] def makePollQuery(workflowId: WorkflowId, replyTo: ActorRef, run: Run) = {
    JesStatusPollQuery(workflowId, replyTo, run)
  }

  private def handleQueryFailure(failure: JesApiQueryFailed) = if (failure.query.failedAttempts < maxRetries) {
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
        val workBatch = JesPollingWorkBatch(work)
        workPullingJesPollingActor ! workBatch
        workInProgress += (workPullingJesPollingActor -> workBatch)
      case None =>
        log.debug("No work for JES poller")
        workPullingJesPollingActor ! NoWorkToDo
    }

    workQueue = beheaded.newWorkQueue
  }

  // Intentionally not final, this runs afoul of SI-4440 (I believe)
  private case class BeheadedWorkQueue(workToDo: Option[NonEmptyList[JesApiQuery]], newWorkQueue: Queue[JesApiQuery])
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
        log.error(s"The JES API worker actor $terminee unexpectedly terminated while conducting ${work.workBatch.tail.size + 1} polls. Making a new one...")
        workInProgress -= terminee
        work.workBatch.toList.foreach {
          case statusQuery: JesStatusPollQuery =>
            self ! JesApiStatusQueryFailed(statusQuery, new JesApiException(throwable))
          case runCreationQuery: JesRunCreationQuery =>
            self ! JesApiRunCreationQueryFailed(runCreationQuery, new JesApiException(throwable))
        }
      case None =>
        // It managed to die while doing absolutely nothing...!?
        // Maybe it deserves an entry in https://en.wikipedia.org/wiki/List_of_unusual_deaths
        // Oh well, in the mean time don't do anything, just start a new one
        log.error(throwable, s"The JES API worker actor managed to unexpectedly terminate whilst doing absolutely nothing (${throwable.getMessage}). This is probably a programming error. Making a new one...")
    }

    resetWorker()
  }

  private def resetWorker() = {
    statusPoller = makeWorkerActor()
    context.watch(statusPoller)
    log.info(s"watching $statusPoller")
  }

  // Separate method to allow overriding in tests:
  private[statuspolling] def makeWorkerActor(): ActorRef = context.actorOf(statusPollerProps)
}

object JesApiQueryManager {

  def props(qps: Int Refined Positive, serviceRegistryActor: ActorRef): Props = Props(new JesApiQueryManager(qps, serviceRegistryActor)).withDispatcher(BackendDispatcher)

  sealed trait JesApiQueryManagerRequest

  /**
    * Poll the job represented by the Run.
    */
  final case class DoPoll(workflowId: WorkflowId, run: Run) extends JesApiQueryManagerRequest

  /**
    * Create an ephemeral pipeline and run it in JES.
    */
  final case class DoCreateRun(workflowId: WorkflowId, genomicsInterface: Genomics, rpr: RunPipelineRequest) extends JesApiQueryManagerRequest

  private[statuspolling] trait JesApiQuery {
    def workflowId: WorkflowId
    val failedAttempts: Int
    def requester: ActorRef
    def genomicsInterface: Genomics
    def withFailedAttempt: JesApiQuery
    def backoff: Backoff
    def httpRequest: HttpRequest
    def contentLength: Long = Option(httpRequest.getContent).map(_.getLength).getOrElse(0L)
  }
  private object JesApiQuery {
    // This must be a def, we want a new one each time (they're mutable! Boo!)
    def backoff: Backoff = SimpleExponentialBackoff(1.second, 1000.seconds, 1.5d)
  }

  private[statuspolling] case class JesStatusPollQuery(workflowId: WorkflowId, requester: ActorRef, run: Run, failedAttempts: Int = 0, backoff: Backoff = JesApiQuery.backoff) extends JesApiQuery {
    override val genomicsInterface = run.genomicsInterface
    override lazy val httpRequest = run.getOperationCommand.buildHttpRequest()
    override def withFailedAttempt = this.copy(failedAttempts = failedAttempts + 1, backoff = backoff.next)
  }

  private[statuspolling] case class JesRunCreationQuery(workflowId: WorkflowId, requester: ActorRef, genomicsInterface: Genomics, rpr: RunPipelineRequest, failedAttempts: Int = 0, backoff: Backoff = JesApiQuery.backoff) extends JesApiQuery {
    override def withFailedAttempt = this.copy(failedAttempts = failedAttempts + 1, backoff = backoff.next)
    override lazy val httpRequest = genomicsInterface.pipelines().run(rpr).buildHttpRequest()
  }

  trait JesApiQueryFailed {
    val query: JesApiQuery
    val cause: JesApiException
  }

  final case class JesApiStatusQueryFailed(query: JesApiQuery, cause: JesApiException) extends JesApiQueryFailed
  final case class JesApiRunCreationQueryFailed(query: JesApiQuery, cause: JesApiException) extends JesApiQueryFailed

  private[statuspolling] final case class JesPollingWorkBatch(workBatch: NonEmptyList[JesApiQuery])
  private[statuspolling] case object NoWorkToDo

  private[statuspolling] final case class RequestJesPollingWork(maxBatchSize: Int)

  final case class GoogleJsonException(e: GoogleJsonError, responseHeaders: HttpHeaders) extends IOException with CromwellFatalExceptionMarker {
    override def getMessage: String = e.getMessage
  }

  class JesApiException(val e: Throwable) extends RuntimeException(e) with CromwellFatalExceptionMarker {
    override def getMessage: String = "Unable to complete JES Api Request"
  }
}
