package cromwell.services.metrics.bard.model
import cromwell.core.WorkflowId

import java.time.Instant

case class TaskSummaryEvent(workflowId: WorkflowId,
                            jobIdKey: String,
                            terminalState: String,
                            cloud: String,
                            dockerImage: String,
                            cpuCount: Int,
                            memoryBytes: Double,
                            startTime: Instant,
                            endTime: Instant
) extends BardEvent {
  override def eventName: String = "task:summary"

}
