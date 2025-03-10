package cromwell.backend.google.batch.api

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import cromwell.backend.google.batch.api.request.{BatchApiRequestHandler, BatchRequestExecutor}
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.services.instrumentation.CromwellInstrumentationActor

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

/**
 * Sends batched requests to Batch as a worker to the BatchApiRequestManager
 */
class BatchApiRequestWorker(val pollingManager: ActorRef,
                            val batchInterval: FiniteDuration,
                            override val serviceRegistryActor: ActorRef,
                            batchRequestExecutor: BatchRequestExecutor
)(implicit val batchHandler: BatchApiRequestHandler)
    extends Actor
    with ActorLogging
    with CromwellInstrumentationActor {

  import BatchApiRequestManager._
  import BatchApiRequestWorker._

  self ! NoWorkToDo // Starts the check-for-work cycle when the actor is fully initialized.

  implicit val ec: ExecutionContext = context.dispatcher

  override def receive = {
    case BatchApiWorkBatch(workBatch) =>
      log.info(s"Got a Batch request batch with ${workBatch.size} requests.")
      handleBatch(workBatch).andThen(interstitialRecombobulation)
      ()
    case NoWorkToDo =>
      scheduleCheckForWork()
  }

  private def handleBatch(workBatch: NonEmptyList[BatchApiRequest]): Future[List[Try[Unit]]] = {
    // Create the batch:
    val batch = workBatch.foldLeft(batchHandler.makeBatchRequest) { case (batch, request) =>
      batchHandler.enqueue(request, batch, pollingManager)
    }

    // Execute the batch and return the list of successes and failures:
    // NB: Blocking and error prone. If this fails, let the supervisor in the BatchApiRequestManager catch and resubmit
    // the work.
    batchRequestExecutor.execute(batch)(ec)
  }

  // TODO: FSMify this actor?
  private def interstitialRecombobulation: PartialFunction[Try[List[Try[Unit]]], Unit] = {
    case Success(allSuccesses) if allSuccesses.forall(_.isSuccess) =>
      log.info(s"All status polls completed successfully.")
      scheduleCheckForWork()
    case Success(someFailures) =>
      val errors = someFailures collect { case Failure(t) => t.getMessage }
      log.warning(
        "Batch request worker had {} failures making {} requests: {}",
        errors.size,
        someFailures.size,
        errors.mkString(System.lineSeparator, "," + System.lineSeparator, "")
      )
      scheduleCheckForWork()
    case Failure(t) =>
      log.info(s"### FIND ME: we have hit programmer error. Stack trace: ${t.getStackTrace.mkString}")
      // NB: Should be impossible since we only ever do completionPromise.trySuccess()
      val msg =
        "Programmer Error in BatchApiRequestWorker.scala: Completion promise unexpectedly set to Failure: {}. Don't do this, otherwise the Future.sequence is short-circuited on the first failure"
      log.error(msg, t.getMessage)
      scheduleCheckForWork()
  }

  /**
   * Schedules a check for work.
   * Warning: Only use this from inside a receive method.
   */
  private def scheduleCheckForWork(): Unit = {
    context.system.scheduler.scheduleOnce(batchInterval) {
      pollingManager ! BatchApiRequestManager.BatchWorkerRequestWork(MaxBatchSize)
    }
    ()
  }
}

object BatchApiRequestWorker {
  def props(pollingManager: ActorRef,
            batchInterval: FiniteDuration,
            serviceRegistryActor: ActorRef,
            batchRequestExecutor: BatchRequestExecutor
  )(implicit
    batchHandler: BatchApiRequestHandler
  ): Props =
    Props(new BatchApiRequestWorker(pollingManager, batchInterval, serviceRegistryActor, batchRequestExecutor))
      .withDispatcher(BackendDispatcher)

  // The Batch API limits does not support batched requests, still, this limit seems reasonable (coming from PAPIv2).
  // See: https://cloud.google.com/batch/quotas
  val MaxBatchSize = 100
}
