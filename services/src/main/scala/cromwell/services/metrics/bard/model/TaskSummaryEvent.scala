package cromwell.services.metrics.bard.model

import java.util.UUID
import scala.jdk.CollectionConverters._

case class TaskSummaryEvent(workflowId: UUID,
                            parentWorkflowId: Option[UUID],
                            rootWorkflowId: UUID,
                            jobTag: String,
                            jobFullyQualifiedName: String,
                            jobIndex: Option[Int],
                            jobAttempt: Int,
                            terminalState: String,
                            cloudPlatform: Option[String],
                            dockerImage: Option[String],
                            cpuCount: Int,
                            memoryBytes: Double,
                            startTime: String,
                            cpuStartTime: Option[String],
                            endTime: String,
                            jobSeconds: Long,
                            cpuSeconds: Option[Long]
) extends BardEvent {
  override def eventName: String = "task:summary"

  override def getProperties: java.util.Map[String, Any] =
    (super.assembleScalaProperties.map {
      case (field, None) => (field, null)
      case (field, Some(value)) => (field, value)
      case t => t
    } + ("distinct_id" -> workflowId)).asJava
}
