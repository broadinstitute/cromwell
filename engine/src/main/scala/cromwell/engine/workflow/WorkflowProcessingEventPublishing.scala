package cromwell.engine.workflow

import java.time.OffsetDateTime

import akka.actor.ActorRef
import common.util.VersionUtil
import cromwell.core.WorkflowId
import cromwell.core.WorkflowProcessingEvents.EventKey.{CromwellId, CromwellVersion, Description, Timestamp}
import cromwell.core.WorkflowProcessingEvents.{DescriptionEventValue, ProcessingEventsKey}
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
}
