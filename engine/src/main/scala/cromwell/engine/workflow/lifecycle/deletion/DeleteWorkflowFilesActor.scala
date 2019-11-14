package cromwell.engine.workflow.lifecycle.deletion

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.RootWorkflowId
import cromwell.engine.workflow.lifecycle.deletion.DeleteWorkflowFilesActor.StartWorkflowFilesDeletion
import cromwell.services.metadata.MetadataService.{GetRootAndSubworkflowOutputs, RootAndSubworkflowOutputsLookupFailure, RootAndSubworkflowOutputsLookupResponse}

class DeleteWorkflowFilesActor(rootWorkflowId: RootWorkflowId,
                               serviceRegistryActor: ActorRef) extends Actor with ActorLogging {

  override def receive: Receive = {
    case StartWorkflowFilesDeletion =>
      serviceRegistryActor ! GetRootAndSubworkflowOutputs(rootWorkflowId)
    case o: RootAndSubworkflowOutputsLookupResponse =>
      println(s"received list of outputs!! Outputs: $o")
    case f: RootAndSubworkflowOutputsLookupFailure =>
      println(s"Something went wrong. Error ${f.reason}")
    case other =>
      log.error(s"Programmer Error! The DeleteWorkflowFilesActor for root workflow ${rootWorkflowId.id.toString} " +
      s"received unexpected message! ($sender sent $other})")
  }
}


object DeleteWorkflowFilesActor {

  // Commands
  sealed trait DeleteWorkflowFilesActorCommand
  object StartWorkflowFilesDeletion extends DeleteWorkflowFilesActorCommand

  sealed trait DeleteWorkflowFilesActorState



  def props(rootWorkflowId: RootWorkflowId, serviceRegistryActor: ActorRef): Props = {
    Props(new DeleteWorkflowFilesActor(rootWorkflowId, serviceRegistryActor))
  }
}
