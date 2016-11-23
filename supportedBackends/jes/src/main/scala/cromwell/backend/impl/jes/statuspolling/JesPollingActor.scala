package cromwell.backend.impl.jes.statuspolling

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.batch.json.JsonBatchCallback
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpHeaders
import com.google.api.services.genomics.model.Operation
import cromwell.backend.impl.jes.Run
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager.{JesPollingWorkBatch, JesStatusPollQuery, NoWorkToDo}
import cromwell.backend.impl.jes.statuspolling.JesPollingActor._

import scala.collection.JavaConversions._
import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.util.{Failure, Success, Try}
import scala.concurrent.duration._

/**
  * Polls JES for status. Pipes the results back (so expect either a RunStatus or a akka.actor.Status.Failure).
  */
class JesPollingActor(pollingManager: ActorRef) extends Actor with ActorLogging {

  // We want to query at just under our fixed JES QPS limit of 20 per second. That should hopefully allow some room at the edges
  // for things like new calls, etc.
  val MaxBatchSize = 100
  val BatchInterval = 5.5.seconds
  self ! NoWorkToDo // Starts the check-for-work cycle

  implicit val ec: ExecutionContext = context.dispatcher

  override def receive = {

    case JesPollingWorkBatch(workBatch) =>
      log.debug(s"Got a polling batch with ${workBatch.tail.size + 1} requests.")
      val batchResultFutures = handleBatch(workBatch)
      val overallFuture = Future.sequence(batchResultFutures.toList)
      overallFuture.andThen(interstitialRecombobulation)
      ()
    case NoWorkToDo =>
      scheduleCheckForWork()
  }

  private def handleBatch(workBatch: NonEmptyList[JesStatusPollQuery]): NonEmptyList[Future[Try[Unit]]] = {

    // Assume that the auth for the first element is also able to query the remaining Runs
    val batch: BatchRequest = createBatch(workBatch.head.run)

    // Create the batch:
    val batchFutures = workBatch map { pollingRequest =>
      val completionPromise = Promise[Try[Unit]]()
      val resultHandler = batchResultHandler(pollingRequest.requester, completionPromise)
      enqueueStatusPollInBatch(pollingRequest.run, batch, resultHandler)
      completionPromise.future
    }

    // Execute the batch and return the map:
    runBatch(batch)
    batchFutures
  }

  // These are separate functions so that the tests can hook in and replace the JES-side stuff
  private[statuspolling] def createBatch(run: Run): BatchRequest = run.genomicsInterface.batch()
  private[statuspolling] def enqueueStatusPollInBatch(run: Run, batch: BatchRequest, resultHandler: JsonBatchCallback[Operation]) = {
    run.getOperationCommand.queue(batch, resultHandler)
  }
  private[statuspolling] def runBatch(batch: BatchRequest) = batch.execute()

  private def batchResultHandler(originalRequester: ActorRef, completionPromise: Promise[Try[Unit]]) = new JsonBatchCallback[Operation] {
    override def onSuccess(operation: Operation, responseHeaders: HttpHeaders): Unit = {
      log.debug(s"Batch result onSuccess callback triggered!")
      originalRequester ! interpretOperationStatus(operation)
      completionPromise.trySuccess(Success(()))
      ()
    }

    override def onFailure(e: GoogleJsonError, responseHeaders: HttpHeaders): Unit = {
      log.debug(s"Batch request onFailure callback triggered!")
      originalRequester ! JesPollFailed(e, responseHeaders)
      completionPromise.trySuccess(Failure(new Exception(mkErrorString(e))))
      ()
    }
  }

  private[statuspolling] def mkErrorString(e: GoogleJsonError) = e.getErrors.toList.mkString(", ")
  private[statuspolling] def interpretOperationStatus(operation: Operation) = Run.interpretOperationStatus(operation)

  // TODO: FSMify this actor?
  private def interstitialRecombobulation: PartialFunction[Try[List[Try[Unit]]], Unit] = {
    case Success(allSuccesses) if allSuccesses.forall(_.isSuccess) =>
      log.debug(s"All status polls completed successfully.")
      scheduleCheckForWork()
    case Success(someFailures) =>
      val errors = someFailures collect { case Failure(t) => t.getMessage }
      if (log.isDebugEnabled) {
        log.warning("{} failures fetching JES statuses", errors.size)
      } else {
        log.warning("{} failures fetching JES statuses: {}", errors.size, errors.mkString(", "))
      }
      scheduleCheckForWork()
    case Failure(t) =>
      // NB: Should be impossible since we only ever do completionPromise.trySuccess()
      log.error("Completion promise unexpectedly set to Failure: {}", t.getMessage)
      scheduleCheckForWork()
  }

  /**
    * Schedules a check for work.
    * Warning: Only use this from inside a receive method.
    */
  private def scheduleCheckForWork(): Unit = {
    context.system.scheduler.scheduleOnce(BatchInterval) { pollingManager ! JesApiQueryManager.RequestJesPollingWork(MaxBatchSize) }
    ()
  }
}

object JesPollingActor {
  def props(pollingManager: ActorRef) = Props(new JesPollingActor(pollingManager))

  final case class JesPollFailed(e: GoogleJsonError, responseHeaders: HttpHeaders)
}
