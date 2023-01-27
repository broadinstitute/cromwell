package cromwell.backend.google.pipelines.batch

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.core.WorkflowId

import scala.concurrent.ExecutionContext



object GcpBatchBackendSingletonActor {
  def props(name: String) = Props(new GcpBatchBackendSingletonActor(name))

  case class BatchRequest(workflowId: WorkflowId, projectId: String, region: String, jobName: String, runtimeAttributes: GcpBatchRuntimeAttributes)
  case class BatchGetJob(jobId: String)
}

final class GcpBatchBackendSingletonActor (name: String) extends Actor with ActorLogging {

  import GcpBatchBackendSingletonActor._

  implicit val ec: ExecutionContext = context.dispatcher

  override def receive: Receive = {
    case jobSubmission: BatchRequest =>
      //val job = GcpBatchJob(jobSubmission, 2000, 200, "e2-standard-4", "gcr.io/google-containers/busybox")
      val job = GcpBatchJob(jobSubmission,200,200, "e2-standard-4", jobSubmission.runtimeAttributes)
      job.submitJob()
      //result.getStatus
    case jobStatus: BatchGetJob =>
      println("matched job status")
      val gcpBatchPoll = new GcpBatchJobGetRequest
      gcpBatchPoll.GetJob(jobStatus.jobId)
      ()
    case other =>
      log.error("Unknown message to GCP Batch Singleton Actor: {}. Dropping it.", other)

  }

}


