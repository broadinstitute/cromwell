package cromwell.backend.impl.jes.statuspolling

import java.io.IOException

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.backend.impl.jes.statuspolling.JesPollingActor.JesPollFailed
import cromwell.backend.impl.jes.{Run, RunStatus}

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success}

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
      pollingActorClientPromise foreach { p =>
        log.debug(s"Polled status received: $r")
        p.complete(Success(r))
      }
      pollingActorClientPromise = None

    case JesPollFailed(e, responseHeaders) =>
      // TODO: Explode? Retry? Fail?
      pollingActorClientPromise foreach { p =>
        log.debug("JES poll failed! Sad.")
        p.complete(Failure(new IOException(s"Google request failed: ${e.toPrettyString}")))
      }
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
