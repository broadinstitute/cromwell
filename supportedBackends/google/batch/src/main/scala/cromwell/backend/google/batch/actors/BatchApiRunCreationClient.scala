package cromwell.backend.google.batch.actors

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.backend.google.batch.api.{BatchApiRequestManager, GcpBatchRequestFactory}
import cromwell.backend.google.batch.api.BatchApiRequestManager.{
  BatchApiRunCreationQueryFailed,
  SystemBatchApiException
}
import cromwell.backend.google.batch.models.GcpBatchRequest
import cromwell.backend.google.batch.monitoring.BatchInstrumentation
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.logging.JobLogger

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

/**
  * Handles the flow for submitting a single job to GCP
  */
trait BatchApiRunCreationClient { this: Actor with ActorLogging with BatchInstrumentation =>
  private var runCreationClientPromise: Option[Promise[StandardAsyncJob]] = None

  // handles messages produced from GcpBatchBackendSingletonActor
  def runCreationClientReceive: Actor.Receive = {
    case job: StandardAsyncJob =>
      log.info(s"runCreationClientReceive -> StandardAsyncJob: ${job.jobId}")
      runSuccess()
      completePromise(Success(job))
    case BatchApiRunCreationQueryFailed(_, e) =>
      log.error(e, s"runCreationClientReceive -> BatchApiRunCreationQueryFailed: ${e.getMessage}")
      completePromise(Failure(e))
  }

  private def completePromise(job: Try[StandardAsyncJob]): Unit = {
    runCreationClientPromise.foreach {
      _.complete(job)
    }
    runCreationClientPromise = None
  }

  def runBatchJob(
    request: GcpBatchRequest,
    backendSingletonActor: ActorRef,
    requestFactory: GcpBatchRequestFactory,
    jobLogger: JobLogger
  ): Future[StandardAsyncJob] =
    runCreationClientPromise match {
      case Some(p) =>
        p.future
      case None =>
        jobLogger.info(s"Asking singleton actor to submit a job: ${request.jobName}")

        backendSingletonActor ! BatchApiRequestManager.BatchRunCreationRequest(
          request.workflowId,
          self,
          requestFactory.submitRequest(request)
        )
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
        new Exception("The job was removed from the queue before a Batch creation request was made")
      )
}
