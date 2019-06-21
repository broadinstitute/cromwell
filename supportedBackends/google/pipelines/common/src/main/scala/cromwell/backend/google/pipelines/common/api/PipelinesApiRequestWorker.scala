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
import scala.util.control.NoStackTrace
import scala.util.{Failure, Success, Try}

/**
  * Sends batched requests to JES as a worker to the JesApiQueryManager
  */
class PipelinesApiRequestWorker(val pollingManager: ActorRef, val batchInterval: FiniteDuration, override val serviceRegistryActor: ActorRef)
                               (implicit val batchHandler: PipelinesApiRequestHandler)
  extends Actor with ActorLogging with CromwellInstrumentationActor {

  // The interval to delay between submitting each batch
  log.info("PAPI request worker batch interval is {}", batchInterval)

  self ! NoWorkToDo // Starts the check-for-work cycle when the actor is fully initialized.

  implicit val ec: ExecutionContext = context.dispatcher

  // Batches can be re-used, and since this actor only executes one at a time, we only need one
  private lazy val batch = batchHandler.makeBatchRequest

  override def receive = {
    case PipelinesApiWorkBatch(workBatch) =>
      log.debug(s"Got a PAPI request batch with ${workBatch.tail.size + 1} requests.")
      handleBatch(workBatch).andThen(interstitialRecombobulation)
      ()
    case NoWorkToDo =>
      scheduleCheckForWork()
  }

  private def handleBatch(workBatch: NonEmptyList[PAPIApiRequest]): Future[List[Try[Unit]]] = {
    // Assume that the auth for the first element is also good enough for everything else:
    val batch: BatchRequest = createBatch()

    // Create the batch:
    // TODO: WARNING: These call change 'batch' as a side effect. Things might go awry if map runs items in parallel?
    val batchFutures = workBatch map { batchHandler.enqueue(_, batch, pollingManager) }

    // Execute the batch and return the list of successes and failures:
    // NB: Blocking and error prone. If this fails, let the supervisor in the PipelinesApiRequestManager catch and resubmit
    // the work.
    runBatch(batch)
    Future.sequence(batchFutures.toList)
  }
  
  // These are separate functions so that the tests can hook in and replace the JES-side stuff
  private[api] def createBatch(): BatchRequest = batch
  private[api] def runBatch(batch: BatchRequest): Unit = {
    try {
      if (batch.size() > 0) batch.execute()
    } catch {
      case e: java.io.IOException =>
        val msg = s"A batch of PAPI status requests failed. The request manager will retry automatically up to 10 times. The error was: ${e.getMessage}"
        throw new Exception(msg, e.getCause) with NoStackTrace
    }
  }

  // TODO: FSMify this actor?
  private def interstitialRecombobulation: PartialFunction[Try[List[Try[Unit]]], Unit] = {
    case Success(allSuccesses) if allSuccesses.forall(_.isSuccess) =>
      log.debug(s"All status polls completed successfully.")
      scheduleCheckForWork()
    case Success(someFailures) =>
      val errors = someFailures collect { case Failure(t) => t.getMessage }
      log.warning("PAPI request worker had {} failures making {} requests: {}", errors.size, someFailures.size, errors.mkString(System.lineSeparator, "," + System.lineSeparator, ""))
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
