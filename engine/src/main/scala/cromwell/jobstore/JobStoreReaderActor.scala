package cromwell.jobstore

import akka.actor.{ActorLogging, ActorRef, Props}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.actor.BatchActor.CommandAndReplyTo
import cromwell.core.actor.ThrottlerActor
import cromwell.jobstore.JobStoreActor.{JobComplete, JobNotComplete, JobStoreReadFailure, QueryJobCompletion}

import scala.util.{Failure, Success}

object JobStoreReaderActor {
  def props(database: JobStore) = Props(new JobStoreReaderActor(database)).withDispatcher(EngineDispatcher)
}

class JobStoreReaderActor(database: JobStore) extends ThrottlerActor[CommandAndReplyTo[QueryJobCompletion]] with ActorLogging {
  override def processHead(head: CommandAndReplyTo[QueryJobCompletion]) = {
    val action = database.readJobResult(head.command.jobKey, head.command.taskOutputs) 
    action onComplete {
      case Success(Some(result)) => head.replyTo ! JobComplete(result)
      case Success(None) => head.replyTo ! JobNotComplete
      case Failure(t) =>
        log.error(t, "JobStoreReadFailure")
        head.replyTo ! JobStoreReadFailure(t)
    }
    action
  }

  override def commandToData(snd: ActorRef) = {
    case query: QueryJobCompletion => CommandAndReplyTo(query, sender())
  }
}
