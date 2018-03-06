package cromwell.backend.impl.jes.statuspolling

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager.JesApiAbortQueryFailed
import cromwell.backend.impl.jes.statuspolling.RunAbort.JesAbortRequestSuccessful
import cromwell.core.logging.JobLogging

trait JesAbortClient { this: Actor with ActorLogging with JobLogging =>
  val pollingActor: ActorRef

  def abortActorClientReceive: Actor.Receive = {
    case JesAbortRequestSuccessful(jobId) =>
      jobLogger.info(s"Successfully requested cancellation of $jobId")
    case JesApiAbortQueryFailed(jobId, e) =>
      jobLogger.error(s"Could not request cancellation of job $jobId", e)
  }
}
