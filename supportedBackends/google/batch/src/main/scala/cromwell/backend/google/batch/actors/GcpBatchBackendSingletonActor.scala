package cromwell.backend.google.batch.actors

import akka.actor.{Actor, ActorLogging, ActorRef, Props, Timers}
import com.google.cloud.batch.v1.JobName
import com.google.longrunning.Operation
import cromwell.backend.BackendSingletonActorAbortWorkflow
import cromwell.backend.google.batch.api.{GcpBatchApiRequestHandler, GcpBatchRequestFactory}
import cromwell.backend.google.batch.models.GcpBatchRequest
import cromwell.backend.google.batch.monitoring.BatchInstrumentation
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.services.instrumentation.CromwellInstrumentationScheduler

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

object GcpBatchBackendSingletonActor {
  def props(requestFactory: GcpBatchRequestFactory, serviceRegistryActor: ActorRef)(implicit requestHandler: GcpBatchApiRequestHandler): Props = {
    Props(new GcpBatchBackendSingletonActor(requestFactory, serviceRegistryActor = serviceRegistryActor))
      .withDispatcher(BackendDispatcher)
  }

  // This is the only type of messages that can be processed by this actor from this actor
  sealed trait Action extends Product with Serializable
  object Action {
    final case class SubmitJob(request: GcpBatchRequest) extends Action
    final case class QueryJob(jobName: JobName) extends Action
    final case class AbortJob(jobName: JobName) extends Action
  }

  // This is the only type of messages produced from this actor while reacting to received messages
  sealed trait Event extends Product with Serializable
  object Event {
    final case class JobSubmitted(job: com.google.cloud.batch.v1.Job) extends Event
    final case class JobStatusRetrieved(job: com.google.cloud.batch.v1.Job) extends Event
    final case class JobAbortRequestSent(operation: Operation) extends Event
    final case class ActionFailed(jobName: String, cause: Throwable) extends Event
  }

}

final class GcpBatchBackendSingletonActor(requestFactory: GcpBatchRequestFactory, override val serviceRegistryActor: ActorRef)(implicit requestHandler: GcpBatchApiRequestHandler)
  extends Actor
    with ActorLogging
    with BatchInstrumentation
    with CromwellInstrumentationScheduler
    with Timers {

  import GcpBatchBackendSingletonActor._

  private implicit val ec: ExecutionContext = context.dispatcher

  override def preStart() = {
    startInstrumentationTimer()
    super.preStart()
  }

  private def normalReceive: Receive = {
    case Action.SubmitJob(request) =>
      val replyTo = sender()
      log.info(s"Submitting job (${request.jobName}) to GCP, workflowId = ${request.workflowId}")
      Future {
        requestHandler.submit(requestFactory.submitRequest(request))
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
          log.info(s"Job ($jobName) retrieved from GCP, state = ${job.getStatus.getState}")
          replyTo ! Event.JobStatusRetrieved(job)

        case Failure(exception) =>
          log.error(exception,  s"Failed to query job status ($jobName) from GCP")
          replyTo ! Event.ActionFailed(jobName.toString ,exception)
      }

    case Action.AbortJob(jobName) =>
      val replyTo = sender()

      Future {
        requestHandler.abort(requestFactory.abortRequest(jobName))
      }.onComplete {
        case Success(operation) =>
          log.info(s"Job ($jobName) aborted from GCP")
          replyTo ! Event.JobAbortRequestSent(operation)

        case Failure(exception) =>
          log.error(exception, s"Failed to abort job ($jobName) from GCP")
          replyTo ! Event.ActionFailed(jobName.toString, exception)
      }

    // Cromwell sends this message
    case BackendSingletonActorAbortWorkflow(workflowId) =>
      // It seems that AbortJob(jobName) is processed before this message, hence, we don't need to do anything else.
      // If it ever becomes necessary, we'll need to create link submitted jobs to its workflow id, which require
      // us to be cautious because batch deletes jobs instead of canceling them, hence, we should not delete jobs
      // that are on a final state.
      log.info(s"Cromwell requested to abort workflow $workflowId")

    case other => log.error(s"Unexpected message from {} to ${this.getClass.getSimpleName}: {}", sender().path.name, other)
  }

  override def receive = instrumentationReceive(loadMetricHandler _).orElse(normalReceive)

  private def loadMetricHandler() = {
    // TODO: Implement this once we have details to expose
//    val load = if (workQueue.size > LoadConfig.PAPIThreshold) HighLoad else NormalLoad
//    serviceRegistryActor ! LoadMetric("PAPIQueryManager", load)
//    updateQueueSize(workQueue.size)
  }
}
