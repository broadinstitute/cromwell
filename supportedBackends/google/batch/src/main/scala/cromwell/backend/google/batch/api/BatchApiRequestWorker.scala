package cromwell.backend.google.batch.api

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import cromwell.backend.google.batch.api.request.{BatchApiRequestHandler, GcpBatchGroupedRequests}
//import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.http.javanet.NetHttpTransport
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.services.instrumentation.CromwellInstrumentationActor

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.control.NoStackTrace
import scala.util.{Failure, Success, Try}

// TODO: Alex - porting this file is almost done
/**
 * Sends batched requests to JES as a worker to the JesApiQueryManager
 */
class BatchApiRequestWorker(val pollingManager: ActorRef,
                            val batchInterval: FiniteDuration,
                            override val serviceRegistryActor: ActorRef
)(implicit val batchHandler: BatchApiRequestHandler)
    extends Actor
    with ActorLogging
    with CromwellInstrumentationActor {

  import BatchApiRequestWorker._
  import BatchApiRequestManager._

  self ! NoWorkToDo // Starts the check-for-work cycle when the actor is fully initialized.

  implicit val ec: ExecutionContext = context.dispatcher

  override def receive = {
    case BatchApiWorkBatch(workBatch) =>
      log.debug(s"Got a Batch request batch with ${workBatch.tail.size + 1} requests.")
      handleBatch(workBatch).andThen(interstitialRecombobulation)
      ()
    case NoWorkToDo =>
      scheduleCheckForWork()
  }

  private def handleBatch(workBatch: NonEmptyList[BatchApiRequest]): Future[List[Try[Unit]]] = {
    // Create the batch:
    val batch = workBatch.foldLeft(createBatch()) { case (batch, request) =>
      batchHandler.enqueue(request, batch, pollingManager)
    }

    // Execute the batch and return the list of successes and failures:
    // NB: Blocking and error prone. If this fails, let the supervisor in the BatchApiRequestManager catch and resubmit
    // the work.
    //
    // TODO: Alex - verify that the work is actually resubmitted on failures
    runBatch(batch)
  }

  // These are separate functions so that the tests can hook in and replace the JES-side stuff
  private[api] def createBatch(): GcpBatchGroupedRequests = batchHandler.makeBatchRequest
  private[api] def runBatch(batch: GcpBatchGroupedRequests): Future[List[Try[Unit]]] = {
    val result =
      if (batch.size > 0) batch.execute(ec)
      else Future.successful(List.empty)

    result.recover {
      // TODO: Alex - this seems unnecessary
      case e: java.io.IOException =>
        val msg =
          s"A batch of Batch status requests failed. The request manager will retry automatically up to 10 times. The error was: ${e.getMessage}"
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
      log.warning(
        "Batch request worker had {} failures making {} requests: {}",
        errors.size,
        someFailures.size,
        errors.mkString(System.lineSeparator, "," + System.lineSeparator, "")
      )
      scheduleCheckForWork()
    case Failure(t) =>
      // NB: Should be impossible since we only ever do completionPromise.trySuccess()
      val msg =
        "Programmer Error: Completion promise unexpectedly set to Failure: {}. Don't do this, otherwise the Future.sequence is short-circuited on the first failure"
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
  def props(pollingManager: ActorRef, batchInterval: FiniteDuration, serviceRegistryActor: ActorRef)(implicit
    batchHandler: BatchApiRequestHandler
  ): Props =
    Props(new BatchApiRequestWorker(pollingManager, batchInterval, serviceRegistryActor))
      .withDispatcher(BackendDispatcher)

  // TODO: Verify this fact applies to GCP Batch
  // https://cloud.google.com/batch/quotas
  // The Batch API limits us to 100 at a time
  val MaxBatchSize = 100

  // TODO: Alex - seems unused
  val HttpTransport: NetHttpTransport = GoogleNetHttpTransport.newTrustedTransport
}
