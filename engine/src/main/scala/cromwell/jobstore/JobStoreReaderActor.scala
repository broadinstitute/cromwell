package cromwell.jobstore

import akka.actor.{Actor, ActorLogging, Props}
import akka.event.LoggingReceive
import cromwell.jobstore.JobStoreActor.{JobComplete, JobNotComplete, JobStoreReadFailure, QueryJobCompletion}

import scala.util.{Failure, Success}

object JobStoreReaderActor {
  def props(database: JobStoreDatabase) = Props(new JobStoreReaderActor(database))
}

class JobStoreReaderActor(database: JobStoreDatabase) extends Actor with ActorLogging {

  implicit val ec = context.dispatcher

  override def receive = LoggingReceive {
    case QueryJobCompletion(key) =>
      val replyTo = sender()
      database.readJobResult(key) onComplete {
        case Success(Some(result)) => replyTo ! JobComplete(result)
        case Success(None) => replyTo ! JobNotComplete
        case Failure(t) => replyTo ! JobStoreReadFailure(t)
      }
    case unknownMessage => log.error(s"Unexpected message to JobStoreReader: $unknownMessage")
  }
}
