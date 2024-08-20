package cromwell.backend.google.pipelines.common.api.clients

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.backend.google.pipelines.common.PapiInstrumentation
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.PipelinesApiStatusQueryFailed
import cromwell.backend.google.pipelines.common.api.{PipelinesApiRequestFactory, PipelinesApiRequestManager, RunStatus}
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.WorkflowId

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

/**
  * I'm putting this stuff in a mixin to avoid polluting the main class.
  *
  * Be sure to make the main class's receive look like:
  *   override def receive = pollingActorClientReceive orElse { ... }
  */
trait PipelinesApiStatusRequestClient { this: Actor with ActorLogging with PapiInstrumentation =>

  private var pollingActorClientPromise: Option[Promise[RunStatus]] = None

  val papiApiActor: ActorRef
  val requestFactory: PipelinesApiRequestFactory

  def pollingActorClientReceive: Actor.Receive = {
    case r: RunStatus =>
      log.debug(s"Polled status received: $r")
      pollSuccess()
      completePromise(Success(r))
    case PipelinesApiStatusQueryFailed(_, e) =>
      log.debug("JES poll failed!")
      completePromise(Failure(e))
  }

  private def completePromise(runStatus: Try[RunStatus]) = {
    pollingActorClientPromise foreach { _.complete(runStatus) }
    pollingActorClientPromise = None
  }

  def pollStatus(workflowId: WorkflowId, jobId: StandardAsyncJob): Future[RunStatus] =
    pollingActorClientPromise match {
      case Some(p) => p.future
      case None =>
        papiApiActor ! PipelinesApiRequestManager.PAPIStatusPollRequest(workflowId,
                                                                        self,
                                                                        requestFactory.getRequest(jobId),
                                                                        jobId
        )
        val newPromise = Promise[RunStatus]()
        pollingActorClientPromise = Option(newPromise)
        newPromise.future
    }
}
