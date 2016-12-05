package cromwell.jobstore

import akka.actor.{Actor, ActorLogging, Props}
import akka.event.LoggingReceive
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.jobstore.JobStoreActor.{JobComplete, JobNotComplete, JobStoreReadFailure, QueryJobCompletion}

import scala.util.{Failure, Success}

object JobStoreReaderActor {
  def props(database: JobStore) = Props(new JobStoreReaderActor(database)).withDispatcher(EngineDispatcher)
}

class JobStoreReaderActor(database: JobStore) extends Actor with ActorLogging {

  implicit val ec = context.dispatcher

  override def receive = LoggingReceive {
    case QueryJobCompletion(key, taskOutputs) =>
      val replyTo = sender()
      database.readJobResult(key, taskOutputs) onComplete {
        case Success(Some(result)) => replyTo ! JobComplete(result)
        case Success(None) => replyTo ! JobNotComplete
        case Failure(t) =>
          log.error(t, "JobStoreReadFailure")
          replyTo ! JobStoreReadFailure(t)
      }
    case unknownMessage => log.error(s"Unexpected message to JobStoreReader: $unknownMessage")
  }
}
