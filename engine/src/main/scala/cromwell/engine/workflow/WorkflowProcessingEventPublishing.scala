package cromwell.engine.workflow

import java.time.OffsetDateTime

import akka.actor.ActorRef
import cats.Monad
import common.util.VersionUtil
import common.validation.IOChecked.IOChecked
import cromwell.core.WorkflowProcessingEvents.EventKey.{CromwellId, CromwellVersion, Description, Timestamp}
import cromwell.core.WorkflowProcessingEvents.{DescriptionEventValue, ProcessingEventsKey}
import cromwell.core.{WorkflowId, WorkflowMetadataKeys}
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}

import scala.util.Random

object WorkflowProcessingEventPublishing {
  private lazy val cromwellVersion = VersionUtil.getVersion("cromwell")

  def publish(workflowId: WorkflowId, cromwellId: String, descriptionValue: DescriptionEventValue.Value, serviceRegistry: ActorRef): Unit = {
    def randomNumberString: String = Random.nextInt(Int.MaxValue).toString

    def metadataKey(workflowId: WorkflowId, randomNumberString: String, key: String) =
      MetadataKey(workflowId = workflowId, jobKey = None, s"$ProcessingEventsKey[$randomNumberString]:$key")

    val random = randomNumberString

    val processingFields = List(
      Description.key -> descriptionValue.value,
      CromwellId.key -> cromwellId,
      Timestamp.key -> OffsetDateTime.now(),
      CromwellVersion.key -> cromwellVersion
    )

    val metadata = processingFields map { case (k, v) =>
      MetadataEvent(metadataKey(workflowId = workflowId, randomNumberString = random, key = k), MetadataValue(v))
    }

    serviceRegistry ! PutMetadataAction(metadata)
  }

  def publishLabelsToMetadata(workflowId: WorkflowId,
                              labels: Map[String, String],
                              serviceRegistry: ActorRef): IOChecked[Unit] = {
    val defaultLabel = "cromwell-workflow-id" -> s"cromwell-$workflowId"
    Monad[IOChecked].pure(labelsToMetadata(workflowId, labels + defaultLabel, serviceRegistry))
  }

  private def labelsToMetadata(workflowId: WorkflowId,
                               labels: Map[String, String],
                               serviceRegistry: ActorRef): Unit = {
    labels foreach { case (labelKey, labelValue) =>
      val metadataKey = MetadataKey(workflowId, None, s"${WorkflowMetadataKeys.Labels}:$labelKey")
      val metadataValue = MetadataValue(labelValue)
      serviceRegistry ! PutMetadataAction(MetadataEvent(metadataKey, metadataValue))
    }
  }
}
