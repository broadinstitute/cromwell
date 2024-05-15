package cromwell.services.metrics.bard.model

import java.util.UUID

case class TaskSummaryEvent(workflowId: UUID,
                            jobIdKey: String,
                            terminalState: String,
                            cloud: String,
                            dockerImage: String,
                            cpuCount: Int,
                            memoryBytes: Double,
                            startTime: String,
                            endTime: String
) extends BardEvent {
  override def eventName: String = "task:summary"

}
