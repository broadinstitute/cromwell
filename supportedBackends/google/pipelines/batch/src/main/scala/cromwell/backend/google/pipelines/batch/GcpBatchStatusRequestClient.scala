package cromwell.backend.google.pipelines.batch

import akka.actor.{Actor, ActorLogging, ActorRef}
import com.google.cloud.batch.v1.{BatchServiceClient, JobName}
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.PipelinesApiStatusQueryFailed
//import cromwell.backend.google.pipelines.common.api.{PipelinesApiRequestManager, RunStatus}
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.WorkflowId

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

trait GcpBatchStatusRequestClient { this: Actor with ActorLogging =>

  private var pollingActorClientPromise: Option[Promise[GcpBatchRunStatus]] = None

  val gcpBatchActor: ActorRef
  //val requestFactory: PipelinesApiRequestFactory

  def pollingActorClientReceive: Actor.Receive = {
    case r: GcpBatchRunStatus =>
      log.debug(s"Polled status received: $r")
      //pollSuccess()
      completePromise(Success(r))
    case PipelinesApiStatusQueryFailed(_, e) => // update for batch
      log.debug("JES poll failed!")
      completePromise(Failure(e))
  }

  private def completePromise(runStatus: Try[GcpBatchRunStatus]) = {
    pollingActorClientPromise foreach { _.complete(runStatus) }
    pollingActorClientPromise = None
  }

  def pollStatus(workflowId: WorkflowId, jobId: StandardAsyncJob): Future[GcpBatchRunStatus] = {
    pollingActorClientPromise match {
      case Some(p) => p.future
      case None =>
        //gcpBatchActor ! PipelinesApiRequestManager.PAPIStatusPollRequest(workflowId, self, requestFactory.getRequest(jobId), jobId)
        gcpBatchActor ! GcpBatchSingleton
        val newPromise = Promise[GcpBatchRunStatus]()
        pollingActorClientPromise = Option(newPromise)
        newPromise.future
    }
  }


}

object GcpBatchSingleton {

  val projectId = "batch-testing-350715"
  val region = "us-central1"
  val jobName = "hello-dan-04"

  val batchServiceClient = BatchServiceClient
    .create()

  val job = batchServiceClient
    .getJob(JobName
      .newBuilder()
      .setProject(projectId)
      .setLocation(region)
      .setJob(jobName)
      .build())

}
