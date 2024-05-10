package cromwell.backend.google.batch.actors

import akka.actor.{Actor, ActorLogging, ActorRef}
import com.google.cloud.batch.v1.JobName
import cromwell.backend.google.batch.api.BatchApiRequestManager.{BatchApiStatusQueryFailed, BatchStatusPollRequest}
import cromwell.backend.google.batch.api.GcpBatchRequestFactory
import cromwell.backend.google.batch.models.RunStatus
import cromwell.backend.google.batch.monitoring.BatchInstrumentation
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.WorkflowId

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

/**
  * Allows fetching a job status
  */
trait BatchApiStatusRequestClient { this: Actor with ActorLogging with BatchInstrumentation =>

  private var pollingActorClientPromise: Option[Promise[RunStatus]] = None

  def pollingActorClientReceive: Actor.Receive = {
    case status: RunStatus =>
      log.debug(s"Polled status received: $status")
      pollSuccess()
      completePromise(Success(status))
    case BatchApiStatusQueryFailed(query, e) =>
      log.error(e, s"Poll status failed for job ${query.jobId}: ${e.getMessage}")
      completePromise(Failure(e))
  }

  private def completePromise(result: Try[RunStatus]): Unit = {
    pollingActorClientPromise foreach { _.complete(result) }
    pollingActorClientPromise = None
  }

  def pollStatus(
    workflowId: WorkflowId,
    jobName: JobName,
    backendSingletonActor: ActorRef,
    requestFactory: GcpBatchRequestFactory
  ): Future[RunStatus] =
    pollingActorClientPromise match {
      case Some(p) => p.future
      case None =>
        backendSingletonActor ! BatchStatusPollRequest(
          workflowId,
          self,
          requestFactory.queryRequest(jobName),
          StandardAsyncJob(jobName.toString)
        )

        val newPromise = Promise[RunStatus]()
        pollingActorClientPromise = Option(newPromise)
        newPromise.future
    }
}
