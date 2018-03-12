package cromwell.backend.impl.jes.statuspolling

import akka.actor.{Actor, ActorLogging, ActorRef}
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model.RunPipelineRequest
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager.{PAPIApiException, JesApiRunCreationQueryFailed}
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.WorkflowId

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

object JesRunCreationClient {

  /**
    * Exception used to represent the fact that a job was aborted before a creation attempt was made.
    * Meaning it was in the queue when the abort request was made, so it was just removed from the queue.
    */
  case object JobAbortedException extends PAPIApiException(new Exception("The job was removed from the queue before a PAPI creation request was made"))
}

/**
  * I'm putting this stuff in a mixin to avoid polluting the main class.
  *
  * Be sure to make the main class's receive look like:
  *   override def receive = runCreationClientReceive orElse { ... }
  */
trait JesRunCreationClient { this: Actor with ActorLogging =>

  private var runCreationClientPromise: Option[Promise[StandardAsyncJob]] = None

  val pollingActor: ActorRef

  def runCreationClientReceive: Actor.Receive = {
    case job: StandardAsyncJob =>
      completePromise(Success(job))
    case JesApiRunCreationQueryFailed(_, e) => completePromise(Failure(e))
  }

  private def completePromise(job: Try[StandardAsyncJob]) = {
    runCreationClientPromise foreach { _.complete(job) }
    runCreationClientPromise = None
  }

  def runPipeline(workflowId: WorkflowId, genomicsInterface: Genomics, rpr: RunPipelineRequest): Future[StandardAsyncJob] = {
    runCreationClientPromise match {
      case Some(p) =>
        p.future
      case None =>
        pollingActor ! JesApiQueryManager.DoCreateRun(workflowId, genomicsInterface, rpr)
        val newPromise = Promise[StandardAsyncJob]()
        runCreationClientPromise = Option(newPromise)
        newPromise.future
    }
  }
}
