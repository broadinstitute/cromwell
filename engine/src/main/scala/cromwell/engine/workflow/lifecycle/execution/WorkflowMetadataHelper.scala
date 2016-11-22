package cromwell.engine.workflow.lifecycle.execution

import java.time.OffsetDateTime

import akka.actor.ActorRef
import cromwell.core.{WorkflowId, WorkflowMetadataKeys, WorkflowState}
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}

import scala.util.Random

trait WorkflowMetadataHelper {

  def serviceRegistryActor: ActorRef
  
  def pushWorkflowStart(workflowId: WorkflowId) = {
    val startEvent = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.StartTime), MetadataValue(OffsetDateTime.now.toString))
    serviceRegistryActor ! PutMetadataAction(startEvent)
  }
  
  def pushWorkflowEnd(workflowId: WorkflowId) = {
    val metadataEventMsg = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.EndTime), MetadataValue(OffsetDateTime.now.toString))
    serviceRegistryActor ! PutMetadataAction(metadataEventMsg)
  }
  
  def pushWorkflowFailures(workflowId: WorkflowId, failures: List[Throwable]) = {
    val failureEvents = failures flatMap { r => throwableToMetadataEvents(MetadataKey(workflowId, None, s"${WorkflowMetadataKeys.Failures}[${Random.nextInt(Int.MaxValue)}]"), r) }
    serviceRegistryActor ! PutMetadataAction(failureEvents)
  }
  
  def pushCurrentStateToMetadataService(workflowId: WorkflowId, workflowState: WorkflowState): Unit = {
    val metadataEventMsg = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.Status),
      MetadataValue(workflowState))
    serviceRegistryActor ! PutMetadataAction(metadataEventMsg)
  }
  
}
