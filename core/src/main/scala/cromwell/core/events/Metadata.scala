package cromwell.core.events

import cromwell.core.WorkflowId

sealed trait MetadataAlert {
  def workflowId: WorkflowId
  def count: Long
}
final case class HeavyMetadataAlert(workflowId: WorkflowId, count: Long) extends MetadataAlert
final case class MaxMetadataAlert(workflowId: WorkflowId, count: Long, limit: Long) extends MetadataAlert
