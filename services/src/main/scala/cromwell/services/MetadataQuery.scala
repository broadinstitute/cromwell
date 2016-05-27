package cromwell.services

import java.sql.Timestamp
import cromwell.core.{KnowsWhatTimeItIs, WorkflowId}
import org.slf4j.LoggerFactory
import wdl4s.values.{WdlBoolean, WdlFloat, WdlInteger, WdlValue}

case class MetadataJobKey(callFqn: String, index: Option[Int], attempt: Int)

case class MetadataKey(workflowId: WorkflowId, jobKey: Option[MetadataJobKey], key: String)

object MetadataEvent extends KnowsWhatTimeItIs {
  def apply(key: MetadataKey, value: MetadataValue) = {
    // TODO: Look into db/slick using OffsetDateTime, or storing datetimes as UTC?
    // http://stackoverflow.com/questions/34608650/scala-slick-3-0-implicit-mapping-between-java8-offsetdatetime-and-timestamp
    // https://github.com/slick/slick/issues/1026
    new MetadataEvent(key, value, currentTime)
  }
}


sealed trait MetadataType { def typeName: String }
case object MetadataString extends MetadataType { override val typeName = "string" }
case object MetadataInt extends MetadataType { override val typeName = "int" }
case object MetadataNumber extends MetadataType { override val typeName = "number" }
case object MetadataBoolean extends MetadataType { override val typeName = "boolean" }

object MetadataValue {
  def apply(value: Any) = {
    value match {
      case null => new MetadataValue(null, null)
      case WdlInteger(i) => new MetadataValue(i.toString, MetadataInt)
      case WdlFloat(f) => new MetadataValue(f.toString, MetadataNumber)
      case WdlBoolean(b) => new MetadataValue(b.toString, MetadataBoolean)
      case value: WdlValue => new MetadataValue(value.valueString, MetadataString)
      case _: Int | Long => new MetadataValue(value.toString, MetadataInt)
      case _: Double | Float => new MetadataValue(value.toString, MetadataNumber)
      case _: Boolean => new MetadataValue(value.toString, MetadataBoolean)
      case _ => new MetadataValue(value.toString, MetadataString)
    }
  }
}

object MetadataType {
  val log = LoggerFactory.getLogger("Metadata Type")

  def fromString(s: String) = s match {
    case MetadataString.typeName => MetadataString
    case MetadataInt.typeName => MetadataInt
    case MetadataNumber.typeName => MetadataNumber
    case MetadataBoolean.typeName => MetadataBoolean
    case _ =>
      log.warn(s"Unknown Metadata type $s. Falling back to MetadataString type")
      MetadataString
  }
}

case class MetadataValue(value: String, valueType: MetadataType)

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
