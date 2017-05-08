package cromwell.backend.impl.jes.statuspolling

import java.io.IOException

import akka.actor.{Actor, ActorLogging, ActorRef, Props, Terminated}
import cats.data.NonEmptyList
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpHeaders
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model.RunPipelineRequest
import cromwell.backend.impl.jes.Run
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager._
import cromwell.core.CromwellFatalExceptionMarker
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core.retry.{Backoff, SimpleExponentialBackoff}
import cromwell.util.StopAndLogSupervisor
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric._

import scala.concurrent.duration._
import scala.collection.immutable.Queue

/**
  * Holds a set of JES API requests until a JesQueryActor pulls the work.
  */
class JesApiQueryManager(val qps: Int Refined Positive) extends Actor with ActorLogging with StopAndLogSupervisor {

  private implicit val ec = context.dispatcher
  private val maxRetries = 10

  // workQueue is protected for the unit tests, not intended to be generally overridden
  protected[statuspolling] var workQueue: Queue[JesApiQuery] = Queue.empty
  private var workInProgress: Map[ActorRef, JesPollingWorkBatch] = Map.empty

  private def statusPollerProps = JesPollingActor.props(self, qps)

  // statusPoller is protected for the unit tests, not intended to be generally overridden
  protected[statuspolling] var statusPoller: ActorRef = _

  resetWorker()

  override def receive = {
    case DoPoll(run) => workQueue :+= JesStatusPollQuery(sender, run)
    case DoCreateRun(genomics, rpr) =>
      workQueue :+= JesRunCreationQuery(sender, genomics, rpr)
    case q: JesApiQuery => workQueue :+= q
    case RequestJesPollingWork(maxBatchSize) =>
      log.debug("Request for JES Polling Work received (max batch: {}, current queue size is {})", maxBatchSize, workQueue.size)
      handleJesPollingRequest(sender, maxBatchSize)
    case failure: JesApiQueryFailed => handleQueryFailure(failure)
    case Terminated(actorRef) => handleTerminated(actorRef)
    case other => log.error(s"Unexpected message to JesPollingManager: $other")
  }

  private def handleQueryFailure(failure: JesApiQueryFailed) = if (failure.query.failedAttempts < maxRetries) {
    val nextRequest = failure.query.withFailedAttempt
    val delay = nextRequest.backoff.backoffMillis.millis
    context.system.scheduler.scheduleOnce(delay, self, nextRequest)
    ()
  } else {
    failure.query.requester ! failure
    ()
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

  private final case class BeheadedWorkQueue(workToDo: Option[NonEmptyList[JesApiQuery]], newWorkQueue: Queue[JesApiQuery])
  private def beheadWorkQueue(maxBatchSize: Int): BeheadedWorkQueue = {

    val head = workQueue.take(maxBatchSize).toList
    val tail = workQueue.drop(maxBatchSize)

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
  }
  private object JesApiQuery {
    // This must be a def, we want a new one each time (they're mutable! Boo!)
    def backoff: Backoff = SimpleExponentialBackoff(1.second, 1000.seconds, 1.5d)
  }
  private[statuspolling] final case class JesStatusPollQuery(requester: ActorRef, run: Run, failedAttempts: Int = 0, backoff: Backoff = JesApiQuery.backoff) extends JesApiQuery {
    override val genomicsInterface = run.genomicsInterface
    override def withFailedAttempt = this.copy(failedAttempts = failedAttempts + 1, backoff = backoff.next)
  }
  private[statuspolling] final case class JesRunCreationQuery(requester: ActorRef, genomicsInterface: Genomics, rpr: RunPipelineRequest, failedAttempts: Int = 0, backoff: Backoff = JesApiQuery.backoff) extends JesApiQuery {
    override def withFailedAttempt = this.copy(failedAttempts = failedAttempts + 1, backoff = backoff.next)
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
