package cromwell.services.metadata

import java.nio.file.Path
import java.time.OffsetDateTime

import cats.data.NonEmptyList
import cromwell.core.WorkflowId
import cromwell.core.path.PathImplicits._
import org.slf4j.LoggerFactory
import wdl4s.values.{WdlBoolean, WdlFloat, WdlInteger, WdlOptionalValue, WdlValue}

case class MetadataJobKey(callFqn: String, index: Option[Int], attempt: Int)

case class MetadataKey(workflowId: WorkflowId, jobKey: Option[MetadataJobKey], key: String)

object MetadataKey {

  val KeySeparator = ':'

  def apply(workflowId: WorkflowId, jobKey: Option[MetadataJobKey], keys: String*): MetadataKey = {
    new MetadataKey(workflowId, jobKey, compositeKey(keys:_*))
  }

  def compositeKey(keys: String*) = keys.toList.mkString(KeySeparator.toString)
}

object MetadataEvent {
  def apply(key: MetadataKey, value: MetadataValue) = new MetadataEvent(key, Option(value), OffsetDateTime.now)
  def empty(key: MetadataKey) = new MetadataEvent(key, None, OffsetDateTime.now)
}


sealed trait MetadataType { def typeName: String }
case object MetadataString extends MetadataType { override val typeName = "string" }
case object MetadataInt extends MetadataType { override val typeName = "int" }
case object MetadataNumber extends MetadataType { override val typeName = "number" }
case object MetadataBoolean extends MetadataType { override val typeName = "boolean" }

object MetadataValue {
  def apply(value: Any): MetadataValue = {
    Option(value).getOrElse("") match {
      case WdlInteger(i) => new MetadataValue(i.toString, MetadataInt)
      case WdlFloat(f) => new MetadataValue(f.toString, MetadataNumber)
      case WdlBoolean(b) => new MetadataValue(b.toString, MetadataBoolean)
      case WdlOptionalValue(_, Some(o)) => apply(o)
      case value: WdlValue => new MetadataValue(value.valueString, MetadataString)
      case _: Int | Long => new MetadataValue(value.toString, MetadataInt)
      case _: Double | Float => new MetadataValue(value.toString, MetadataNumber)
      case _: Boolean => new MetadataValue(value.toString, MetadataBoolean)
      case path: Path => new MetadataValue(path.toRealString, MetadataString)
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

case class MetadataEvent(key: MetadataKey, value: Option[MetadataValue], offsetDateTime: OffsetDateTime)

case class MetadataQueryJobKey(callFqn: String, index: Option[Int], attempt: Int)

object MetadataQueryJobKey {
  def forMetadataJobKey(jobKey: MetadataJobKey) = MetadataQueryJobKey(jobKey.callFqn, jobKey.index, jobKey.attempt)
}

case class MetadataQuery(workflowId: WorkflowId, jobKey: Option[MetadataQueryJobKey], key: Option[String],
                         includeKeysOption: Option[NonEmptyList[String]],
                         excludeKeysOption: Option[NonEmptyList[String]],
                         expandSubWorkflows: Boolean)

object MetadataQuery {
  def forWorkflow(workflowId: WorkflowId) = MetadataQuery(workflowId, None, None, None, None, expandSubWorkflows = false)

  def forJob(workflowId: WorkflowId, jobKey: MetadataJobKey) = {
    MetadataQuery(workflowId, Option(MetadataQueryJobKey.forMetadataJobKey(jobKey)), None, None, None, expandSubWorkflows = false)
  }

  def forKey(key: MetadataKey) = {
    MetadataQuery(key.workflowId, key.jobKey map MetadataQueryJobKey.forMetadataJobKey, Option(key.key), None, None, expandSubWorkflows = false)
  }
}
