package cromwell.backend.google.pipelines.batch

import akka.actor.{Actor, ActorLogging, Props}
import scala.concurrent.ExecutionContext



object GcpBatchBackendSingletonActor {
  def props(name: String) = Props(new GcpBatchBackendSingletonActor(name))

  case class BatchRequest(projectId: String, region: String, jobName: String)

}

final class GcpBatchBackendSingletonActor (name: String) extends Actor with ActorLogging {

  import GcpBatchBackendSingletonActor._

  implicit val ec: ExecutionContext = context.dispatcher


  def receive: Receive = {
    case jobSubmission: BatchRequest =>
      val job = GCPBatchJob(jobSubmission, 2000, 200, "e2-standard-4", "gcr.io/google-containers/busybox")
      job.submitJob()
    case other =>
      log.error("Unknown message to GCP Batch Singleton Actor: {}. Dropping it.", other)

  }

}


