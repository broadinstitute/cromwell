package cromwell.backend.google.batch.actors

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.backend.google.batch.api.BatchApiRequestManager.SystemBatchApiException
import cromwell.backend.google.batch.models.GcpBatchRequest
import cromwell.backend.google.batch.monitoring.BatchInstrumentation
import cromwell.backend.standard.StandardAsyncJob

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

/**
  * Handles the flow for submitting a single job to GCP, we can't do anything when that fails
  */
trait BatchApiRunCreationClient { this: Actor with ActorLogging with BatchInstrumentation =>
  private var runCreationClientPromise: Option[Promise[StandardAsyncJob]] = None

  // handles messages produced from GcpBatchBackendSingletonActor
  def runCreationClientReceive: Actor.Receive = {
    case GcpBatchBackendSingletonActor.Event.JobSubmitted(job) =>
      log.info(s"Job submitted to GCP: ${job.getName}")
      runSuccess()
      completePromise(Success(StandardAsyncJob(job.getName)))

    case GcpBatchBackendSingletonActor.Event.ActionFailed(jobName, cause) =>
      val msg = s"Failed to submit job ($jobName) to GCP"
      log.error(cause, msg)
      runFailed()
      completePromise(Failure(cause))
  }

  private def completePromise(job: Try[StandardAsyncJob]): Unit = {
    runCreationClientPromise.foreach {
      _.complete(job)
    }
    runCreationClientPromise = None
  }

  def runBatchJob(request: GcpBatchRequest, backendSingletonActor: ActorRef): Future[StandardAsyncJob] =
    runCreationClientPromise match {
      case Some(p) =>
        p.future
      case None =>
        log.info(s"Asking singleton actor to submit a job: ${request.jobName}")
        backendSingletonActor ! GcpBatchBackendSingletonActor.Action.SubmitJob(request)
        val newPromise = Promise[StandardAsyncJob]()
        runCreationClientPromise = Option(newPromise)
        newPromise.future
    }
}

object BatchApiRunCreationClient {

  /**
   * Exception used to represent the fact that a job was aborted before a creation attempt was made.
   * Meaning it was in the queue when the abort request was made, so it was just removed from the queue.
   */
  case object JobAbortedException
      extends SystemBatchApiException(
        new Exception("The job was removed from the queue before a PAPI creation request was made")
      )
}
