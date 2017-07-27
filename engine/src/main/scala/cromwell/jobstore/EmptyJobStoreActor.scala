package cromwell.jobstore

import akka.actor.{Actor, Props}
import cromwell.jobstore.JobStoreActor._
import cromwell.core.Dispatcher.EngineDispatcher

class EmptyJobStoreActor extends Actor {
  override def receive: Receive = {
    case w: JobStoreWriterCommand => sender ! JobStoreWriteSuccess(w)
    case _: QueryJobCompletion => sender ! JobNotComplete
  }
}

object EmptyJobStoreActor {
  def props: Props = Props(new EmptyJobStoreActor()).withDispatcher(EngineDispatcher)
}
