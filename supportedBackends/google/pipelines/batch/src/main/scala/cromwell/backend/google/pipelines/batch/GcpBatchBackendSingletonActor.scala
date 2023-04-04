package cromwell.backend.google.pipelines.batch

import akka.actor.{Actor, ActorLogging, Props}
import com.google.cloud.batch.v1.{BatchServiceClient, GetJobRequest, JobName}
import cromwell.core.Dispatcher.BackendDispatcher

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

object GcpBatchBackendSingletonActor {
  def props(name: String) = Props(new GcpBatchBackendSingletonActor(name)).withDispatcher(BackendDispatcher)

  // This is the only type of messages that can be processed by this actor from this actor
  sealed trait Action extends Product with Serializable
  object Action {
    final case class SubmitJob(request: GcpBatchRequest) extends Action
    final case class QueryJobStatus(jobId: String, projectId: String, region: String) extends Action
  }

  // This is the only type of messages produced from this actor while reacting to received messages
  sealed trait Event extends Product with Serializable
  object Event {
    final case class JobSubmitted(job: com.google.cloud.batch.v1.Job) extends Event
    final case class JobStatusRetrieved(job: com.google.cloud.batch.v1.Job) extends Event
    final case class ActionFailed(cause: Throwable) extends Event
  }

}

// TODO: Alex - serviceRegistryActor required to do instrumentation
final class GcpBatchBackendSingletonActor (name: String) extends Actor with ActorLogging {

  import GcpBatchBackendSingletonActor._

  implicit val ec: ExecutionContext = context.dispatcher

  // TODO: Alex - Find the usefulness for the actor because we are not holding any mutable state in it
  override def receive: Receive = {
    case Action.SubmitJob(request) =>
      // TODO: Alex - Consider not hardcoding machineType
      val job = GcpBatchJob(request,"n1-standard-4")
      val replyTo = sender()
      log.info(s"Submitting job (${request.jobName}) to GCP, workflowId = ${request.workflowId}")
      Future {
        job.submitJob()
      }.onComplete {
        case Failure(exception) =>
          log.error(exception, s"Failed to submit job (${request.jobName}) to GCP, workflowId = ${request.workflowId}")
          replyTo ! Event.ActionFailed(exception)

        case Success(job) =>
          log.info(s"Job (${request.jobName}) submitted to GCP, workflowId = ${request.workflowId}, id = ${job.getUid}")
          replyTo ! Event.JobSubmitted(job)
      }

    case Action.QueryJobStatus(jobId, projectId, region) =>
      val replyTo = sender()
      val jobName = JobName.of(projectId, region, jobId)

      // TODO: Alex - Consider creating this client only once in the app lifecycle
      val batchServiceClient = BatchServiceClient.create
      val resultF = Future {
        val request = GetJobRequest.newBuilder.setName(jobName.toString).build
        batchServiceClient.getJob(request)
      }

      // close client
      resultF.onComplete(_ => batchServiceClient.close())
      resultF.onComplete {
        case Success(job) =>
          log.info(s"Job status ($jobName) retrieved from GCP, state = ${job.getStatus.getState}")
          replyTo ! Event.JobStatusRetrieved(job)

        case Failure(exception) =>
          log.error(exception,  s"Failed to query job status ($jobName) from GCP")
          replyTo ! Event.ActionFailed(exception)
      }

    case other =>
      log.error("Unknown message to GCP Batch Singleton Actor: {}. Dropping it.", other)

  }

}