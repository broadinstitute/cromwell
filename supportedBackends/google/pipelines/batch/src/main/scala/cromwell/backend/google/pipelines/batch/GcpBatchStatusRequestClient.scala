package cromwell.backend.google.pipelines.batch

import akka.actor.{Actor, ActorLogging, ActorRef}
//import cromwell.backend.google.pipelines.batch.GcpBatchBackendSingletonActor.GcpBatchJobSuccess
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.PipelinesApiStatusQueryFailed
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.WorkflowId

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

trait GcpBatchStatusRequestClient { this: Actor with ActorLogging =>

  private var pollingActorClientPromise: Option[Promise[RunStatus]] = None

  val gcpBatchActor: ActorRef
  //val requestFactory: PipelinesApiRequestFactory



  def pollingActorClientReceive: Actor.Receive = {
    case r: RunStatus =>
      log.debug(s"Polling status received: $r")
      //pollSuccess()
      completePromise(Success(r))
    case PipelinesApiStatusQueryFailed(_, e) => // update for batch
      log.debug("GCP Batch poll failed!")
      completePromise(Failure(e))
    //case j: GcpBatchJobSuccess =>
        //log.info("received job success")
        //completePromise(Success(j))
    //case other => println(f"poll receive test $other")
  }

  private def completePromise(runStatus: Try[RunStatus]) = {
    pollingActorClientPromise foreach { _.complete(runStatus) }
    pollingActorClientPromise = None
  }

  def pollStatus(workflowId: WorkflowId, jobId: StandardAsyncJob, gcpBatchJobId: String): Future[RunStatus] = {

    val test = new GcpBatchJobGetRequest

    pollingActorClientPromise match {
      case Some(p) =>
        p.future
      case None =>
        println("polling with gcp status request client")
        //gcpBatchActor ! GcpBatchSingleton
        gcpBatchActor ! test.GetJob(gcpBatchJobId)

        println(RunStatus)
        val newPromise = Promise[RunStatus]()
        pollingActorClientPromise = Option(newPromise)
        newPromise.future
    }
  }

}
