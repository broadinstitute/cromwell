package cromwell.backend.google.pipelines.batch

import akka.actor.{Actor, ActorLogging, Props}
//import akka.stream.impl.StreamSubscriptionTimeoutSupport.CancelingSubscriber.onComplete
//import akka.stream.impl.VirtualProcessor.Inert.subscriber.onComplete
//import com.google.cloud.batch.v1.JobStatus
import cromwell.core.WorkflowId

//import scala.concurrent.duration._
import scala.concurrent.ExecutionContext
//import scala.concurrent.{ExecutionContext, Promise}



object GcpBatchBackendSingletonActor {
  def props(name: String) = Props(new GcpBatchBackendSingletonActor(name))

  case class BatchRequest(workflowId: WorkflowId, projectId: String, region: String, jobName: String, runtimeAttributes: GcpBatchRuntimeAttributes)
  case class BatchGetJob(jobId: String)
  case class BatchJobAsk(test: String)
}

final class GcpBatchBackendSingletonActor (name: String) extends Actor with ActorLogging {

  import GcpBatchBackendSingletonActor._

  implicit val ec: ExecutionContext = context.dispatcher

  override def receive: Receive = {
    case jobSubmission: BatchRequest =>
      //val job = GcpBatchJob(jobSubmission, 2000, 200, "e2-standard-4", "gcr.io/google-containers/busybox")
      val job = GcpBatchJob(jobSubmission,200,200, "e2-standard-4", jobSubmission.runtimeAttributes)
      job.submitJob()
      sender() ! job
      //result.getStatus
    case jobStatus: BatchGetJob =>
      println("matched job status")
      println(jobStatus.jobId)
      val gcpBatchPoll = new GcpBatchJobGetRequest
      gcpBatchPoll.GetJob(jobStatus.jobId)
      ()
    case _: BatchJobAsk =>
      println("matched job ask")
      sender() ! "ask"
    case other =>
      log.error("Unknown message to GCP Batch Singleton Actor: {}. Dropping it.", other)

  }

}


