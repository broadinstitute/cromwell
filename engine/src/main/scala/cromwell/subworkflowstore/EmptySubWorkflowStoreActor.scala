package cromwell.subworkflowstore

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.subworkflowstore.SubWorkflowStoreActor._
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

class EmptySubWorkflowStoreActor extends Actor with ActorLogging {
  override def receive: Receive = {
    case register: RegisterSubWorkflow => sender() ! SubWorkflowStoreRegisterSuccess(register)
    case query: QuerySubWorkflow => sender() ! SubWorkflowNotFound(query)
    case _: WorkflowComplete => // No-op!
    case ShutdownCommand => context stop self
    case unknown => log.error(s"SubWorkflowStoreActor received unknown message: $unknown")
  }
}

object EmptySubWorkflowStoreActor {
  def props: Props = Props(new EmptySubWorkflowStoreActor()).withDispatcher(EngineDispatcher)
}
