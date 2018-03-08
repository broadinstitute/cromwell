package cromwell.backend.impl.jes.statuspolling

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.services.genomics.Genomics
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager._
import cromwell.backend.impl.jes.statuspolling.JesPollingActor._
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.services.instrumentation.CromwellInstrumentationActor

import scala.collection.JavaConverters._
import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

/**
  * Sends batched requests to JES as a worker to the JesApiQueryManager
  */
class JesPollingActor(val pollingManager: ActorRef, val batchInterval: FiniteDuration, override val serviceRegistryActor: ActorRef) extends Actor with ActorLogging
  with StatusPolling with RunCreation with RunAbort with CromwellInstrumentationActor {
  // The interval to delay between submitting each batch
  log.info("JES batch polling interval is {}", batchInterval)

  self ! NoWorkToDo // Starts the check-for-work cycle when the actor is fully initialized.

  implicit val ec: ExecutionContext = context.dispatcher

  override def receive = {

    case JesPollingWorkBatch(workBatch) =>
      log.debug(s"Got a polling batch with ${workBatch.tail.size + 1} requests.")
      handleBatch(workBatch).andThen(interstitialRecombobulation)
      ()
    case NoWorkToDo =>
      scheduleCheckForWork()
  }

  private def handleBatch(workBatch: NonEmptyList[PAPIApiRequest]): Future[List[Try[Unit]]] = {
    // Assume that the auth for the first element is also good enough for everything else:
    val batch: BatchRequest = createBatch(workBatch.head.genomicsInterface)

    // Create the batch:
    // WARNING: These call change 'batch' as a side effect. Things might go awry if map runs items in parallel?
    val batchFutures = workBatch map {
      case pollingRequest: PAPIStatusPollRequest => enqueueStatusPollInBatch(pollingRequest, batch)
      case runCreationRequest: PAPIRunCreationRequest => enqueueRunCreationInBatch(runCreationRequest, batch)
      case abortRequest: PAPIAbortRequest => enqueueAbortInBatch(abortRequest, batch)

      // We do the "successful Failure" thing so that the Future.sequence doesn't short-out immediately when the first one fails.
      case other => Future.successful(Failure(new RuntimeException(s"Cannot handle ${other.getClass.getSimpleName} requests")))
    }

    // Execute the batch and return the map:
    runBatch(batch)
    Future.sequence(batchFutures.toList)
  }

  // These are separate functions so that the tests can hook in and replace the JES-side stuff
  private[statuspolling] def createBatch(genomicsInterface: Genomics): BatchRequest = genomicsInterface.batch()
  private[statuspolling] def runBatch(batch: BatchRequest) = batch.execute()

  // TODO: FSMify this actor?
  private def interstitialRecombobulation: PartialFunction[Try[List[Try[Unit]]], Unit] = {
    case Success(allSuccesses) if allSuccesses.forall(_.isSuccess) =>
      log.debug(s"All status polls completed successfully.")
      scheduleCheckForWork()
    case Success(someFailures) =>
      val errors = someFailures collect { case Failure(t) => t.getMessage }
      if (log.isDebugEnabled) {
        log.warning("{} failures (from {} requests) fetching JES statuses", errors.size, someFailures.size)
      } else {
        log.warning("{} failures (from {} requests) fetching JES statuses: {}", errors.size, someFailures.size, errors.mkString(", "))
      }
      scheduleCheckForWork()
    case Failure(t) =>
      // NB: Should be impossible since we only ever do completionPromise.trySuccess()
      val msg = "Programmer Error: Completion promise unexpectedly set to Failure: {}. Don't do this, otherwise the Future.sequence is short-circuited on the first failure"
      log.error(msg, t.getMessage)
      scheduleCheckForWork()
  }

  /**
    * Schedules a check for work.
    * Warning: Only use this from inside a receive method.
    */
  private def scheduleCheckForWork(): Unit = {
    context.system.scheduler.scheduleOnce(batchInterval) { pollingManager ! JesApiQueryManager.RequestJesPollingWork(MaxBatchSize) }
    ()
  }

  private[statuspolling] def mkErrorString(e: GoogleJsonError) = e.getErrors.asScala.toList.mkString(", ")
}

object JesPollingActor {
  def props(pollingManager: ActorRef, batchInterval: FiniteDuration, serviceRegistryActor: ActorRef) = {
    Props(new JesPollingActor(pollingManager, batchInterval, serviceRegistryActor)).withDispatcher(BackendDispatcher)
  }

  // The Batch API limits us to 100 at a time
  val MaxBatchSize = 100
}
