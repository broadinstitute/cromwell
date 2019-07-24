package cromwell.services.metadata.impl.hybridcarbonite

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.core.WorkflowId
import cromwell.services.metadata.MetadataService.{GetLabels, GetLogs, GetMetadataQueryAction, GetSingleWorkflowMetadataAction, GetStatus, MetadataWriteAction, PutMetadataAction, PutMetadataActionAndRespond, ReadAction, WorkflowOutputs, WorkflowQuery}

class HybridCarboniteMetadataServiceActor extends Actor with ActorLogging {
  val standardMetadataActor: ActorRef = ???
  val carbonitingMetadataServiceActor: ActorRef = ??? // FIXME: Needs to be defined

  /*
    FIXME: Not a fan of all the copy/paste going on in the receive but we'll need to add a few layers of inheritance in the metadata action hierarchy to clean this up

    FIXME: Also, this receive function is currently blocking everywhere, that will need to be fixed
   */
  override def receive = {
    case message: MetadataWriteAction =>
      // FIXME: Will potentially need to subset the MetadataEvents and forward appropriately

      /*
        FIXME, not sure the best path here, but the following would at least work
           - Partition the events iterable into label updates & non-label updates
           - Create a new MetadataWriteAction of original subtype w/ the non-label updates & send to standard metadata
           - Partition the label updates based on if their underlying workflowid has been carbonited
           - Send a MetadataWriteAction (of original subtype) to standard metadata for non-carbonited
           - Send a MetadataWriteAction (of original subtype) to standard metadata for arbonited
       */

      ???
    case x: GetLabels => standardMetadataActor.forward(x)
    case x: GetStatus => standardMetadataActor.forward(x)
    case x: WorkflowQuery => standardMetadataActor.forward(x)
    case x: GetLogs =>
      if (hasWorkflowBeenCarbonited(x.workflowId)) carbonitingMetadataServiceActor.forward(x)
      else standardMetadataActor.forward(x)
    case x: WorkflowOutputs =>
      if (hasWorkflowBeenCarbonited(x.workflowId)) carbonitingMetadataServiceActor.forward(x)
      else standardMetadataActor.forward(x)
    case x: GetMetadataQueryAction =>
      if (hasWorkflowBeenCarbonited(x.key.workflowId)) carbonitingMetadataServiceActor.forward(x)
      else standardMetadataActor.forward(x)
    case x: GetSingleWorkflowMetadataAction =>
      if (hasWorkflowBeenCarbonited(x.workflowId)) carbonitingMetadataServiceActor.forward(x)
      else standardMetadataActor.forward(x)
    case _ => standardMetadataActor.forward(_)
  }

  private def hasWorkflowBeenCarbonited(workflowId: WorkflowId): Boolean = {
    ???
  }
}
