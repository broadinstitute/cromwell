package cromwell.backend.impl.jes.statuspolling

import java.io.IOException

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.backend.impl.jes.statuspolling.JesPollingActor.JesPollFailed
import cromwell.backend.impl.jes.{Run, RunStatus}

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

/**
  * I'm putting this stuff in a mixin to avoid polluting the main class.
  *
  * Be sure to make the main class's receive look like:
  *   override def receive = pollingActorClientReceive orElse { ... }
  */
trait JesPollingActorClient { this: Actor with ActorLogging =>

  private var pollingActorClientPromise: Option[Promise[RunStatus]] = None

  val pollingActor: ActorRef

  def pollingActorClientReceive: Actor.Receive = {
    case r: RunStatus =>
      log.debug(s"Polled status received: $r")
      completePromise(Success(r))
    case JesPollFailed(e, responseHeaders) =>
      log.debug("JES poll failed! Sad.")
      completePromise(Failure(new IOException(s"Google request failed: ${e.toPrettyString}")))
  }

  private def completePromise(runStatus: Try[RunStatus]) = {
    pollingActorClientPromise foreach { _.complete(runStatus) }
    pollingActorClientPromise = None
  }

  def pollStatus(run: Run): Future[RunStatus] = {
    pollingActorClientPromise match {
      case Some(p) => p.future
      case None =>
        pollingActor ! JesApiQueryManager.DoPoll(run)
        val newPromise = Promise[RunStatus]()
        pollingActorClientPromise = Option(newPromise)
        newPromise.future
    }
  }
}
