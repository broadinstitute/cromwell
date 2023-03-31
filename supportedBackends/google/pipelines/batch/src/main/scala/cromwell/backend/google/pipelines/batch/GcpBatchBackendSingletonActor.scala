package cromwell.backend.google.pipelines.batch

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.core.WorkflowId
import scala.concurrent.ExecutionContext
import cromwell.backend.google.pipelines.batch.GcpBatchRequestFactory.CreatePipelineParameters

object GcpBatchBackendSingletonActor {
  def props(name: String) = Props(new GcpBatchBackendSingletonActor(name))

  case class GcpBatchRequest(workflowId: WorkflowId, createParameters: CreatePipelineParameters,jobName: String, gcpBatchCommand: String, gcpBatchParameters: CreateGcpBatchParameters)
  case class BatchGetJob(jobId: String, projectId: String, region: String)
  case class BatchJobAsk(test: String)

}

final class GcpBatchBackendSingletonActor (name: String) extends Actor with ActorLogging {

  import GcpBatchBackendSingletonActor._

  implicit val ec: ExecutionContext = context.dispatcher

  override def receive: Receive = {
    case jobSubmission: GcpBatchRequest =>
      val job = GcpBatchJob(jobSubmission,"n1-standard-4")
      job.submitJob()
    case jobStatus: BatchGetJob =>
      log.info("matched job status")
      log.info(jobStatus.jobId)
      val gcpBatchPoll = new GcpBatchJobGetRequest
      gcpBatchPoll.GetJob(jobStatus.jobId, jobStatus.projectId, jobStatus.region)
      ()
    case _: BatchJobAsk =>
      log.info("matched job ask")
      sender() ! "Singleton Actor ready!"
    case other =>
      log.error("Unknown message to GCP Batch Singleton Actor: {}. Dropping it.", other)

  }

}