package cromwell.subworkflowstore

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.subworkflowstore.SubWorkflowStoreActor._

class EmptySubWorkflowStoreActor extends Actor with ActorLogging {
  override def receive: Receive = {
    case register: RegisterSubWorkflow => sender() ! SubWorkflowStoreRegisterSuccess(register) 
    case query: QuerySubWorkflow => sender() ! SubWorkflowNotFound(query)
    case complete: WorkflowComplete =>sender() ! SubWorkflowStoreCompleteSuccess(complete)
    case unknown => log.error(s"SubWorkflowStoreActor received unknown message: $unknown")
  }
}

object EmptySubWorkflowStoreActor {
  def props: Props = Props(new EmptySubWorkflowStoreActor())
}
