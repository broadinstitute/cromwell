package cromwell.backend.google.pipelines.common.api

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestWorker._
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.services.instrumentation.CromwellInstrumentationActor

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

/**
  * Sends batched requests to JES as a worker to the JesApiQueryManager
  */
class PipelinesApiRequestWorker(val pollingManager: ActorRef, val batchInterval: FiniteDuration, override val serviceRegistryActor: ActorRef)
                               (implicit val batchHandler: PipelinesApiRequestHandler)
  extends Actor with ActorLogging with CromwellInstrumentationActor {
  // The interval to delay between submitting each batch
  log.info("JES batch polling interval is {}", batchInterval)

  self ! NoWorkToDo // Starts the check-for-work cycle when the actor is fully initialized.

  implicit val ec: ExecutionContext = context.dispatcher

  // Batches can be re-used, and since this actor only executes one at a time, we only need one
  private lazy val batch = batchHandler.makeBatchRequest

  override def receive = {
    case PipelinesApiWorkBatch(workBatch) =>
      log.debug(s"Got a polling batch with ${workBatch.tail.size + 1} requests.")
      handleBatch(workBatch).andThen(interstitialRecombobulation)
      ()
    case NoWorkToDo =>
      scheduleCheckForWork()
  }

  private def handleBatch(workBatch: NonEmptyList[PAPIApiRequest]): Future[List[Try[Unit]]] = {
    // Assume that the auth for the first element is also good enough for everything else:
    val batch: BatchRequest = createBatch()

    // Create the batch:
    // WARNING: These call change 'batch' as a side effect. Things might go awry if map runs items in parallel?
    val batchFutures = workBatch map { batchHandler.enqueue(_, batch, pollingManager) }

    // Execute the batch and return the map:
    runBatch(batch)
    Future.sequence(batchFutures.toList)
  }
  
  // These are separate functions so that the tests can hook in and replace the JES-side stuff
  private[api] def createBatch(): BatchRequest = batch
  private[api] def runBatch(batch: BatchRequest) = if (batch.size() > 0) batch.execute()

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
    context.system.scheduler.scheduleOnce(batchInterval) { pollingManager ! PipelinesApiRequestManager.PipelinesWorkerRequestWork(MaxBatchSize) }
    ()
  }
}

object PipelinesApiRequestWorker {
  def props(pollingManager: ActorRef, batchInterval: FiniteDuration, serviceRegistryActor: ActorRef)
              (implicit batchHandler: PipelinesApiRequestHandler) = {
    Props(new PipelinesApiRequestWorker(pollingManager, batchInterval, serviceRegistryActor)).withDispatcher(BackendDispatcher)
  }

  // The Batch API limits us to 100 at a time
  val MaxBatchSize = 100

  val HttpTransport = GoogleNetHttpTransport.newTrustedTransport
}
