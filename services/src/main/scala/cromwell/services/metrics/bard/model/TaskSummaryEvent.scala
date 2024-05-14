package cromwell.services.metrics.bard.model
import cromwell.core.WorkflowId

import java.util

case class TaskSummaryEvent(workflowId: WorkflowId,
                            jobDescriptorKey: util.Map[String, Any],
                            taskMetadata: util.Map[String, Any]
) extends BardEvent {
  override def eventName: String = "task:summary"

}
