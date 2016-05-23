package cromwell.services

import java.sql.Timestamp
import java.time.OffsetDateTime

import cromwell.core.WorkflowId

case class MetadataJobKey(callFqn: String, index: Option[Int], attempt: Int)

case class MetadataKey(workflowId: WorkflowId, jobKey: Option[MetadataJobKey], key: String)

case class MetadataValue(value: String)

object MetadataEvent {
  def apply(key: MetadataKey, value: MetadataValue) = {
    // TODO: Look into db/slick using OffsetDateTime, or storing datetimes as UTC?
    // http://stackoverflow.com/questions/34608650/scala-slick-3-0-implicit-mapping-between-java8-offsetdatetime-and-timestamp
    // https://github.com/slick/slick/issues/1026
    new MetadataEvent(key, value, Timestamp.valueOf(OffsetDateTime.now.toLocalDateTime))
  }
}

case class MetadataEvent(key: MetadataKey, value: MetadataValue, timestamp: Timestamp)

case class MetadataQueryJobKey(callFqn: String, index: Option[Int], attempt: Int)

object MetadataQueryJobKey {
  def forMetadataJobKey(jobKey: MetadataJobKey) = MetadataQueryJobKey(jobKey.callFqn, jobKey.index, jobKey.attempt)
}

case class MetadataQuery(workflowId: WorkflowId, jobKey: Option[MetadataQueryJobKey], key: Option[String])

object MetadataQuery {
  def forWorkflow(workflowId: WorkflowId) = MetadataQuery(workflowId, None, None)

  def forJob(workflowId: WorkflowId, jobKey: MetadataJobKey) = MetadataQuery(workflowId, Option(MetadataQueryJobKey.forMetadataJobKey(jobKey)), None)

  def forKey(key: MetadataKey) = MetadataQuery(key.workflowId, key.jobKey map MetadataQueryJobKey.forMetadataJobKey, Option(key.key))
}
