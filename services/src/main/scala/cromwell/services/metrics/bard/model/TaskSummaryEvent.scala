package cromwell.services.metrics.bard.model
import cromwell.core.WorkflowId

import java.time.Instant
import java.util

case class TaskSummaryEvent(workflowId: WorkflowId,
                            jobDescriptorKey: util.Map[String, Any],
                            cloud: String,
                            dockerImage: String,
                            cpuCount: Int,
                            memoryBytes: Double,
                            startTime: Instant,
                            endTime: Instant
) extends BardEvent {
  override def eventName: String = "task:summary"

}
