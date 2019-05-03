package cromwell.backend.google.pipelines.common.api.clients

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.backend.google.pipelines.common.PapiInstrumentation
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.{PipelinesApiRunCreationQueryFailed, SystemPAPIApiException}
import cromwell.backend.google.pipelines.common.api.{PipelinesApiRequestFactory, PipelinesApiRequestManager}
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.WorkflowId
import cromwell.core.logging.JobLogger

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

object PipelinesApiRunCreationClient {

  /**
    * Exception used to represent the fact that a job was aborted before a creation attempt was made.
    * Meaning it was in the queue when the abort request was made, so it was just removed from the queue.
    */
  case object JobAbortedException extends SystemPAPIApiException(new Exception("The job was removed from the queue before a PAPI creation request was made"))
}

/**
  * I'm putting this stuff in a mixin to avoid polluting the main class.
  *
  * Be sure to make the main class's receive look like:
  *   override def receive = runCreationClientReceive orElse { ... }
  */
trait PipelinesApiRunCreationClient { this: Actor with ActorLogging with PapiInstrumentation =>
  private var runCreationClientPromise: Option[Promise[StandardAsyncJob]] = None

  val papiApiActor: ActorRef
  val requestFactory: PipelinesApiRequestFactory

  def runCreationClientReceive: Actor.Receive = {
    case job: StandardAsyncJob =>
      runSuccess()
      completePromise(Success(job))
    case PipelinesApiRunCreationQueryFailed(_, e) => completePromise(Failure(e))
  }

  private def completePromise(job: Try[StandardAsyncJob]) = {
    runCreationClientPromise foreach { _.complete(job) }
    runCreationClientPromise = None
  }

  def runPipeline(workflowId: WorkflowId, createPipelineParameters: CreatePipelineParameters, jobLogger: JobLogger): Future[StandardAsyncJob] = {
    runCreationClientPromise match {
      case Some(p) =>
        p.future
      case None =>
        papiApiActor ! PipelinesApiRequestManager.PAPIRunCreationRequest(workflowId, self, requestFactory.runRequest(createPipelineParameters, jobLogger))
        val newPromise = Promise[StandardAsyncJob]()
        runCreationClientPromise = Option(newPromise)
        newPromise.future
    }
  }
}
