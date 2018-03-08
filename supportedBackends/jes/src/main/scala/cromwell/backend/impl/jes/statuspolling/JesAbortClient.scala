package cromwell.backend.impl.jes.statuspolling

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager.JesApiAbortQueryFailed
import cromwell.backend.impl.jes.statuspolling.RunAbort.{PAPIAbortRequestSuccessful, PAPIOperationAlreadyCancelled, PAPIOperationHasAlreadyFinished}
import cromwell.core.logging.JobLogging

trait JesAbortClient { this: Actor with ActorLogging with JobLogging =>
  val pollingActor: ActorRef

  def abortActorClientReceive: Actor.Receive = {
    case PAPIAbortRequestSuccessful(jobId) =>
      jobLogger.info(s"Successfully requested cancellation of $jobId")
    // In this case we could immediately return an aborted handle and spare ourselves a d round of polling
    case PAPIOperationAlreadyCancelled(jobId) =>
      jobLogger.info(s"Operation $jobId was already cancelled")
    // In this case we could immediately return an aborted handle and spare ourselves a d round of polling
    case PAPIOperationHasAlreadyFinished(jobId) =>
      jobLogger.info(s"Operation $jobId has already finished")
    case JesApiAbortQueryFailed(jobId, e) =>
      jobLogger.error(s"Could not request cancellation of job $jobId", e)
  }
}
