package cromwell.backend.google.pipelines.batch

import akka.actor.{Actor, ActorLogging, Props}
import com.google.cloud.batch.v1.JobName
import cromwell.backend.BackendSingletonActorAbortWorkflow
import cromwell.backend.google.pipelines.batch.api.GcpBatchRequestFactory
import cromwell.core.Dispatcher.BackendDispatcher

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

object GcpBatchBackendSingletonActor {
  def props(requestFactory: GcpBatchRequestFactory)(implicit requestHandler: GcpBatchApiRequestHandler) = Props(new GcpBatchBackendSingletonActor(requestFactory))
    .withDispatcher(BackendDispatcher)

  // This is the only type of messages that can be processed by this actor from this actor
  sealed trait Action extends Product with Serializable
  object Action {
    final case class SubmitJob(request: GcpBatchRequest) extends Action
    final case class QueryJob(jobName: JobName) extends Action
  }

  // This is the only type of messages produced from this actor while reacting to received messages
  sealed trait Event extends Product with Serializable
  object Event {
    final case class JobSubmitted(job: com.google.cloud.batch.v1.Job) extends Event
    final case class JobStatusRetrieved(job: com.google.cloud.batch.v1.Job) extends Event
    final case class ActionFailed(jobName: String, cause: Throwable) extends Event
  }

}

// TODO: Alex - serviceRegistryActor required to do instrumentation
final class GcpBatchBackendSingletonActor(requestFactory: GcpBatchRequestFactory)(implicit requestHandler: GcpBatchApiRequestHandler) extends Actor with ActorLogging {

  import GcpBatchBackendSingletonActor._

  implicit val ec: ExecutionContext = context.dispatcher

  // TODO: Alex - Find the usefulness for the actor because we are not holding any mutable state in it
  override def receive: Receive = {
    case Action.SubmitJob(request) =>
      val replyTo = sender()
      log.info(s"Submitting job (${request.jobName}) to GCP, workflowId = ${request.workflowId}")
      Future {
        // TODO: Alex - Consider not hardcoding machineType
        requestHandler.submit(requestFactory.submitRequest("n1-standard-4", request))
      }.onComplete {
        case Failure(exception) =>
          log.error(exception, s"Failed to submit job (${request.jobName}) to GCP, workflowId = ${request.workflowId}")
          replyTo ! Event.ActionFailed(request.jobName, exception)

        case Success(job) =>
          log.info(s"Job (${request.jobName}) submitted to GCP, workflowId = ${request.workflowId}, id = ${job.getUid}")
          replyTo ! Event.JobSubmitted(job)
      }

    case Action.QueryJob(jobName) =>
      val replyTo = sender()

      Future {
        requestHandler.query(requestFactory.queryRequest(jobName))
      }.onComplete {
        case Success(job) =>
          log.info(s"Job status ($jobName) retrieved from GCP, state = ${job.getStatus.getState}")
          replyTo ! Event.JobStatusRetrieved(job)

        case Failure(exception) =>
          log.error(exception,  s"Failed to query job status ($jobName) from GCP")
          replyTo ! Event.ActionFailed(jobName.toString ,exception)
      }

    // Cromwell sends this message
    case BackendSingletonActorAbortWorkflow(workflowId) =>
      // TODO: Alex - how do we handle this?
      log.info(s"Cromwell requested to abort workflow $workflowId")
  }
}
