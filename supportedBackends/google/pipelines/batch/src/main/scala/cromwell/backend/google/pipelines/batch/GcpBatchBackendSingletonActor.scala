package cromwell.backend.google.pipelines.batch

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.core.WorkflowId
import scala.concurrent.ExecutionContext

object GcpBatchBackendSingletonActor {
  def props(name: String) = Props(new GcpBatchBackendSingletonActor(name))

  case class BatchRequest(workflowId: WorkflowId, projectId: String, region: String, jobName: String, runtimeAttributes: GcpBatchRuntimeAttributes, gcpBatchCommand: String)
  case class BatchGetJob(jobId: String)
  case class BatchJobAsk(test: String)

}

final class GcpBatchBackendSingletonActor (name: String) extends Actor with ActorLogging {

  import GcpBatchBackendSingletonActor._

  implicit val ec: ExecutionContext = context.dispatcher

  override def receive: Receive = {
    case jobSubmission: BatchRequest =>
      val job = GcpBatchJob(jobSubmission,200,200, "e2-standard-4", jobSubmission.runtimeAttributes)
      job.submitJob()
    case jobStatus: BatchGetJob =>
      log.info("matched job status")
      log.info(jobStatus.jobId)
      val gcpBatchPoll = new GcpBatchJobGetRequest
      gcpBatchPoll.GetJob(jobStatus.jobId)
      ()
    case _: BatchJobAsk =>
      log.info("matched job ask")
      sender() ! "Singleton Actor ready!"
    case other =>
      log.error("Unknown message to GCP Batch Singleton Actor: {}. Dropping it.", other)

  }

}