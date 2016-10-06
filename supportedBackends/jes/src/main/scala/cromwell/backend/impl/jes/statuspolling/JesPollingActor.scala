package cromwell.backend.impl.jes.statuspolling

import akka.actor.{Actor, ActorLogging, ActorRef, PoisonPill, Props}
import cats.data.NonEmptyList
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.batch.json.JsonBatchCallback
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpHeaders
import com.google.api.services.genomics.model.Operation
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager.{JesPollingWorkBatch, JesStatusPollQuery, NoWorkToDo}
import cromwell.backend.impl.jes.statuspolling.JesPollingActor._

import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.util.{Failure, Success, Try}
import scala.concurrent.duration._

/**
  * Polls JES for status. Pipes the results back (so expect either a RunStatus or a akka.actor.Status.Failure).
  */
class JesPollingActor(pollingManager: ActorRef) extends Actor with ActorLogging {

  implicit val ec: ExecutionContext = context.dispatcher

  override def receive = {

    case JesPollingWorkBatch(workBatch) =>
      log.debug(s"Got a polling batch with ${workBatch.tail.size + 1} requests.")
      val batchResultFutures = handleBatch(workBatch)
      val overallFuture = Future.sequence(batchResultFutures.map(_._2).toList)
      overallFuture.andThen(interstitialRecombobulation)
    case NoWorkToDo =>
      context.system.scheduler.scheduleOnce(10.seconds) { checkForWork() }
  }


  // TODO: safety first! That .queue(...) can throw an IOException!
  private def handleBatch(workBatch: NonEmptyList[JesStatusPollQuery]): NonEmptyList[(ActorRef, Future[Unit])] = {

    // Assume that the first element is also
    val batch: BatchRequest = workBatch.head.run.genomicsInterface.batch()

    // Create the batch:
    val requesterMap = workBatch map { pollingRequest =>
      val completionPromise = Promise[Unit]()
      pollingRequest.run.getOperationCommand.queue(batch, batchResultHandler(pollingRequest, completionPromise))
      (pollingRequest.requester, completionPromise.future)
    }

    // Execute the batch and return the map:
    batch.execute()
    requesterMap
  }

  private def batchResultHandler(originalRequest: JesStatusPollQuery, completionPromise: Promise[Unit]) = new JsonBatchCallback[Operation] {
    override def onSuccess(operation: Operation, responseHeaders: HttpHeaders): Unit = {
      log.debug(s"Batch result onSuccess callback triggered!")
      originalRequest.requester ! originalRequest.run.interpretOperationStatus(operation)
      completionPromise.trySuccess(())
    }

    override def onFailure(e: GoogleJsonError, responseHeaders: HttpHeaders): Unit = {
      log.debug(s"Batch request onFailure callback triggered!")
      originalRequest.requester ! JesPollFailed(e, responseHeaders)
      completionPromise.tryFailure(new Exception()) // TODO: We could recover and retry this?
    }
  }

  // TODO: FSMify this actor!
  private def interstitialRecombobulation: PartialFunction[Try[_], Unit] = {
    case Success(_) =>
      log.debug(s"All status polls completed successfully.")
      checkForWork()
    case Failure(t) =>
      log.error("Error fetching JES status: {}. This JES polling actor is stopping itself")
      self ! PoisonPill
  }

  private def checkForWork() = {
    pollingManager ! JesApiQueryManager.RequestJesPollingWork
  }
}

object JesPollingActor {
  def props(pollingManager: ActorRef) = Props(new JesPollingActor(pollingManager))

  final case class JesPollFailed(e: GoogleJsonError, responseHeaders: HttpHeaders)
}
