package cromwell.backend.google.batch.actors

import akka.actor.{Actor, ActorLogging, ActorRef}
import com.google.cloud.batch.v1.{Job, JobName}
import cromwell.backend.google.batch.api.BatchApiRequestManager.{BatchApiStatusQueryFailed, BatchStatusPollRequest}
import cromwell.backend.google.batch.api.GcpBatchRequestFactory
import cromwell.backend.google.batch.monitoring.BatchInstrumentation
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.WorkflowId

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

/**
  * Allows fetching a job
  */
trait BatchApiFetchJobClient { this: Actor with ActorLogging with BatchInstrumentation =>

  private var pollingActorClientPromise: Option[Promise[Job]] = None

  def pollingActorClientReceive: Actor.Receive = {
    case job: Job =>
      log.debug(s"Polled status received: ${job.getStatus}")
      pollSuccess()
      completePromise(Success(job))
    case BatchApiStatusQueryFailed(_, e) =>
      log.debug("JES poll failed!")
      completePromise(Failure(e))
  }
  // handles messages produced from GcpBatchBackendSingletonActor
//  def pollingActorClientReceive: Actor.Receive = {
//    case GcpBatchBackendSingletonActor.Event.JobStatusRetrieved(job) =>
//      log.info(s"Job retrieved from GCP: ${job.getName}: ${job.getStatus}")
//      pollSuccess()
//      completePromise(Success(job))
//
//    case GcpBatchBackendSingletonActor.Event.ActionFailed(jobName, cause) =>
//      val msg = s"Failed to query job ($jobName) from GCP"
//      log.error(cause, msg)
//      pollFailed()
//      completePromise(Failure(cause))
//  }

  private def completePromise(result: Try[Job]): Unit = {
    pollingActorClientPromise foreach { _.complete(result) }
    pollingActorClientPromise = None
  }

  def fetchJob(workflowId: WorkflowId,
               jobName: JobName,
               backendSingletonActor: ActorRef,
               requestFactory: GcpBatchRequestFactory
  ): Future[Job] =
    pollingActorClientPromise match {
      case Some(p) => p.future
      case None =>
        backendSingletonActor ! BatchStatusPollRequest(
          workflowId,
          self,
          requestFactory.queryRequest(jobName),
          StandardAsyncJob(jobName.toString)
        )
//        backendSingletonActor ! GcpBatchBackendSingletonActor.Action.QueryJob(jobName)

        val newPromise = Promise[Job]()
        pollingActorClientPromise = Option(newPromise)
        newPromise.future
    }
}
