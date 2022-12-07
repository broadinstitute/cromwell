package cromwell.backend.google.pipelines.batch

import akka.actor.{Actor, ActorLogging, ActorRef}
import scala.concurrent.{Future, Promise}
import scala.util.{Success, Try}
import cromwell.backend.standard.StandardAsyncJob
import cromwell.backend.BackendJobDescriptor

trait GcpBatchRunCreationClient {
  // set to none.  Review
  this: Actor with ActorLogging => private var runCreationClientPromise: Option[Promise[StandardAsyncJob]] = None

  val gcpBatchApiActor: ActorRef
  //def gcpBatchApiActor: ActorRef
  //def requestFactory: GcpBatchBackendLifecycleFactory

  val jobDescriptor: BackendJobDescriptor
  lazy val batchJob: GcpBatchJob = GcpBatchJob(jobDescriptor)

  def runCreationClientReceive: Actor.Receive = {
    case job: StandardAsyncJob =>
      //runSuccess() # PAPI. may need to replace
      completePromise(Success(job))
  }

  private def completePromise(job: Try[StandardAsyncJob]) = {
    runCreationClientPromise foreach { _.complete(job) }
    runCreationClientPromise = None
  }

  //def runPipeline(workflowId: WorkflowId, jobLogger: JobLogger): Future[StandardAsyncJob] = {
  def runPipeline(): Future[StandardAsyncJob] = {
      runCreationClientPromise match {
        case Some(p) =>
          p.future
        case None =>
          gcpBatchApiActor ! batchJob.callClient
          val newPromise = Promise[StandardAsyncJob]()
          runCreationClientPromise = Option(newPromise)
          newPromise
            .future
      }
    }


}
