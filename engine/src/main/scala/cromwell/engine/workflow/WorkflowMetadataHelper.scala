package cromwell.engine.workflow

import java.time.OffsetDateTime

import akka.actor.ActorRef
import cromwell.core._
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import cromwell.services.metadata.MetadataService._

trait WorkflowMetadataHelper {

  def serviceRegistryActor: ActorRef
  
  def pushWorkflowStart(workflowId: WorkflowId) = {
    val startEvent = MetadataEvent(
      MetadataKey(workflowId, None, WorkflowMetadataKeys.StartTime),
      MetadataValue(OffsetDateTime.now)
    )
    serviceRegistryActor ! PutMetadataAction(startEvent)
  }
  
  def pushWorkflowEnd(workflowId: WorkflowId) = {
    val metadataEventMsg = MetadataEvent(
      MetadataKey(workflowId, None, WorkflowMetadataKeys.EndTime),
      MetadataValue(OffsetDateTime.now)
    )
    serviceRegistryActor ! PutMetadataAction(metadataEventMsg)
  }
  
  def pushWorkflowFailures(workflowId: WorkflowId, failures: List[Throwable]) = {
    val failureEvents = failures flatMap { r => throwableToMetadataEvents(MetadataKey(workflowId, None, s"${WorkflowMetadataKeys.Failures}"), r) }
    serviceRegistryActor ! PutMetadataAction(failureEvents)
  }
  
  def pushCurrentStateToMetadataService(workflowId: WorkflowId, workflowState: WorkflowState, confirmTo: Option[ActorRef] = None): Unit = {
    val metadataEventMsg = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.Status), MetadataValue(workflowState))

    confirmTo match {
      case None => serviceRegistryActor ! PutMetadataAction(metadataEventMsg)
      case Some(actorRef) => serviceRegistryActor ! PutMetadataActionAndRespond(List(metadataEventMsg), actorRef)
    }
  }
  
}
