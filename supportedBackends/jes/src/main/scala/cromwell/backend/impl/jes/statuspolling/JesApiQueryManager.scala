package cromwell.backend.impl.jes.statuspolling

import java.io.IOException

import akka.actor.{Actor, ActorLogging, ActorRef, Props, Terminated}
import cats.data.NonEmptyList
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.{HttpHeaders, HttpRequest}
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model.RunPipelineRequest
import cromwell.backend.impl.jes.Run
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager.{JesApiException, _}
import cromwell.core.CromwellFatalExceptionMarker
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core.retry.{Backoff, SimpleExponentialBackoff}
import cromwell.util.StopAndLogSupervisor
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric._

import scala.annotation.tailrec
import scala.concurrent.duration._
import scala.collection.immutable.Queue

/**
  * Holds a set of JES API requests until a JesQueryActor pulls the work.
  */
class JesApiQueryManager(val qps: Int Refined Positive) extends Actor with ActorLogging with StopAndLogSupervisor {

  private implicit val ec = context.dispatcher
  private val maxRetries = 10
  // FIXME: Determined empirically that we start to get errors when the batch request approaches 15MB.
  // Find out what the real value is
  private val maxBatchRequestSize = 14 * 1024 * 1024
  private val requestTooLargeException = new JesApiException(new IllegalArgumentException(
    "The task run request has exceeded the maximum PAPI request size." +
      "If you have a task with a very large number of inputs and / or outputs in your workflow you should try to reduce it. " +
      "Depending on your case you could: 1) Zip your input files together and unzip them in the command. 2) Use a file of file names " +
      "and localize the files yourself."
  ))

  // workQueue is protected for the unit tests, not intended to be generally overridden
  protected[statuspolling] var workQueue: Queue[JesApiQuery] = Queue.empty
  private var workInProgress: Map[ActorRef, JesPollingWorkBatch] = Map.empty

  private def statusPollerProps = JesPollingActor.props(self, qps)

  // statusPoller is protected for the unit tests, not intended to be generally overridden
  protected[statuspolling] var statusPoller: ActorRef = _

  resetWorker()

  override def receive = {
    case DoPoll(run) => workQueue :+= makePollQuery(sender, run)
    case DoCreateRun(genomics, rpr) =>
      val creationQuery = makeCreateQuery(sender, genomics, rpr)

      if (creationQuery.contentLength > maxBatchRequestSize) {
        creationQuery.requester ! JesApiRunCreationQueryFailed(creationQuery, requestTooLargeException)
      } else workQueue :+= creationQuery
    case q: JesApiQuery => workQueue :+= q
    case RequestJesPollingWork(maxBatchSize) =>
      log.debug("Request for JES Polling Work received (max batch: {}, current queue size is {})", maxBatchSize, workQueue.size)
      handleJesPollingRequest(sender, maxBatchSize)
    case failure: JesApiQueryFailed => handleQueryFailure(failure)
    case Terminated(actorRef) => handleTerminated(actorRef)
    case other => log.error(s"Unexpected message to JesPollingManager: $other")
  }

  private [statuspolling] def makeCreateQuery(replyTo: ActorRef, genomics: Genomics, rpr: RunPipelineRequest) = {
    JesRunCreationQuery(replyTo, genomics, rpr)
  }

  private [statuspolling] def makePollQuery(replyTo: ActorRef, run: Run) = {
    JesStatusPollQuery(replyTo, run)
  }

  private def handleQueryFailure(failure: JesApiQueryFailed) = if (failure.query.failedAttempts < maxRetries) {
    val nextRequest = failure.query.withFailedAttempt
    val delay = nextRequest.backoff.backoffMillis.millis
    context.system.scheduler.scheduleOnce(delay, self, nextRequest)
    ()
  } else {
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

    /*
      * Keep taking from the head of the queue, making sure it stays under maxBatchSize as well as maxBatchRequestSize.
      * Assumes that each query in the queue can fit in an empty batch (queries with a size > maxBatchSize should no be added to the queue)
      *
      * This is a naive approach where we keep adding queries to the head as long as they don't overflow it, and then we stop.
      * This could leave us with half empty batches and is not very efficient.
      *
      * A somewhat less naive approach would be to keep looking in the queue for another query that would fit.
      *
      * More generally this is the knapsack problem (https://en.wikipedia.org/wiki/Knapsack_problem) which is NP-Complete.
      * If this needs further optimization, finding an approximation algorithm that fits well this case would be a place to start.
      */
    @tailrec
    def behead(queue: Queue[JesApiQuery], head: Vector[JesApiQuery]): Vector[JesApiQuery]  = queue.headOption match {
      case Some(query) if head.size < maxBatchSize && head.map(_.contentLength).sum + query.contentLength < maxBatchRequestSize =>
        behead(queue.tail, head :+ query)
      case _ => head
    }

    // Initialize the head with workQueue.head to avoid deadlocking if the first element is > maxBatchRequestSize, 
    // which should NOT happen as they should not be inserted in the queue in the first place
    val head = if (workQueue.isEmpty) Vector.empty else behead(workQueue.tail, workQueue.headOption.toVector).toList
    val tail = workQueue.drop(head.size)

    head match {
      case h :: t => BeheadedWorkQueue(Option(NonEmptyList(h, t)), tail)
      case Nil => BeheadedWorkQueue(None, Queue.empty)
    }
  }

  private def handleTerminated(terminee: ActorRef) = {
    val cause = getFailureCause(terminee).getOrElse(new RuntimeException("No failure reason recorded"))
    // We assume this is a polling actor. Might change in a future update:
    workInProgress.get(terminee) match {
      case Some(work) =>
        // Most likely due to an unexpected HTTP error, push the work back on the queue and keep going
        log.error(s"The JES API worker actor $terminee unexpectedly terminated while conducting ${work.workBatch.tail.size + 1} polls. Making a new one...")
        workInProgress -= terminee
        work.workBatch.toList.foreach {
          case statusQuery: JesStatusPollQuery =>
            self ! JesApiStatusQueryFailed(statusQuery, new JesApiException(cause))
          case runCreationQuery: JesRunCreationQuery =>
            self ! JesApiRunCreationQueryFailed(runCreationQuery, new JesApiException(cause))
        }
      case None =>
        // It managed to die while doing absolutely nothing...!?
        // Maybe it deserves an entry in https://en.wikipedia.org/wiki/List_of_unusual_deaths
        // Oh well, in the mean time don't do anything, just start a new one
        log.error(cause, s"The JES API worker actor managed to unexpectedly terminate whilst doing absolutely nothing (${cause.getMessage}). This is probably a programming error. Making a new one...")
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

  def props(qps: Int Refined Positive): Props = Props(new JesApiQueryManager(qps)).withDispatcher(BackendDispatcher)

  sealed trait JesApiQueryManagerRequest

  /**
    * Poll the job represented by the Run.
    */
  final case class DoPoll(run: Run) extends JesApiQueryManagerRequest

  /**
    * Create an ephemeral pipeline and run it in JES.
    */
  final case class DoCreateRun(genomicsInterface: Genomics, rpr: RunPipelineRequest) extends JesApiQueryManagerRequest

  private[statuspolling] trait JesApiQuery {
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

  private[statuspolling] case class JesStatusPollQuery(requester: ActorRef, run: Run, failedAttempts: Int = 0, backoff: Backoff = JesApiQuery.backoff) extends JesApiQuery {
    override val genomicsInterface = run.genomicsInterface
    override lazy val httpRequest = run.getOperationCommand.buildHttpRequest()
    override def withFailedAttempt = this.copy(failedAttempts = failedAttempts + 1, backoff = backoff.next)
  }

  private[statuspolling] case class JesRunCreationQuery(requester: ActorRef, genomicsInterface: Genomics, rpr: RunPipelineRequest, failedAttempts: Int = 0, backoff: Backoff = JesApiQuery.backoff) extends JesApiQuery {
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

  final class JesApiException(e: Throwable) extends RuntimeException(e) with CromwellFatalExceptionMarker {
    override def getMessage: String = "Unable to complete JES Api Request"
  }
}
