package cromwell.backend.google.batch.actors

import akka.actor.{Actor, ActorLogging, ActorRef}
import com.google.cloud.batch.v1.{Job, JobName}
import cromwell.backend.google.batch.monitoring.BatchInstrumentation

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

/**
  * Allows fetching a job
  */
trait BatchApiFetchJobClient { this: Actor with ActorLogging with BatchInstrumentation =>

  private var pollingActorClientPromise: Option[Promise[Job]] = None

  // handles messages produced from GcpBatchBackendSingletonActor
  def pollingActorClientReceive: Actor.Receive = {
    case GcpBatchBackendSingletonActor.Event.JobStatusRetrieved(job) =>
      log.info(s"Job retrieved from GCP: ${job.getName}: ${job.getStatus}")
      pollSuccess()
      completePromise(Success(job))

    case GcpBatchBackendSingletonActor.Event.ActionFailed(jobName, cause) =>
      val msg = s"Failed to query job ($jobName) from GCP"
      log.error(cause, msg)
      pollFailed()
      completePromise(Failure(cause))
  }

  private def completePromise(result: Try[Job]): Unit = {
    pollingActorClientPromise foreach { _.complete(result) }
    pollingActorClientPromise = None
  }

  def fetchJob(jobName: JobName, backendSingletonActor: ActorRef): Future[Job] =
    pollingActorClientPromise match {
      case Some(p) => p.future
      case None =>
        backendSingletonActor ! GcpBatchBackendSingletonActor.Action.QueryJob(jobName)

        val newPromise = Promise[Job]()
        pollingActorClientPromise = Option(newPromise)
        newPromise.future
    }
}
