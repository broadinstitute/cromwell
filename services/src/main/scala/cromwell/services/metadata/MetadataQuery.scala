package cromwell.services.metadata

import java.time.OffsetDateTime

import cats.data.NonEmptyList
import cromwell.core._
import cromwell.core.labels.Labels
import org.slf4j.{Logger, LoggerFactory}
import wom.values._
import common.util.TimeUtil._

case class MetadataJobKey(callFqn: String, index: Option[Int], attempt: Int)

case class MetadataKey private (workflowId: WorkflowId, jobKey: Option[MetadataJobKey], key: String)

object MetadataKey {

  val KeySeparator = ':'

  def apply(workflowId: WorkflowId, jobKey: Option[MetadataJobKey], keys: String*): MetadataKey = {
    new MetadataKey(workflowId, jobKey, compositeKey(keys:_*))
  }

  def compositeKey(keys: String*): String = keys.toList.mkString(KeySeparator.toString)
}

object MetadataEvent {
  def apply(key: MetadataKey, value: MetadataValue) = new MetadataEvent(key, Option(value), OffsetDateTime.now)
  def apply(key: MetadataKey, optionalValue: Option[MetadataValue]) = new MetadataEvent(key, optionalValue, OffsetDateTime.now)
  def empty(key: MetadataKey) = new MetadataEvent(key, None, OffsetDateTime.now)

  def labelsToMetadataEvents(labels: Labels, workflowId: WorkflowId): Iterable[MetadataEvent] = {
    labels.value map { label =>
      MetadataEvent(
        MetadataKey(workflowId, None, s"${WorkflowMetadataKeys.Labels}:${label.key}"),
        MetadataValue(label.value)
      )
    }
  }
}

sealed trait MetadataType { def typeName: String }
case object MetadataString extends MetadataType { override val typeName = "string" }
case object MetadataInt extends MetadataType { override val typeName = "int" }
case object MetadataNumber extends MetadataType { override val typeName = "number" }
case object MetadataBoolean extends MetadataType { override val typeName = "boolean" }
/* TODO Might be better to have MetadataNull be a value instead of a type ?
  * We'd need to reorganize MetadataValue, maybe like spray does and have explicit case classes types
  * instead of one generic MetadataValue(value, type)
*/
case object MetadataNull extends MetadataType { override val typeName = "null" }

object MetadataValue {
  def apply(value: Any): MetadataValue = {
    Option(value).getOrElse("") match {
      case WomInteger(i) => new MetadataValue(i.toString, MetadataInt)
      case WomFloat(f) => new MetadataValue(f.toString, MetadataNumber)
      case WomLong(f) => new MetadataValue(f.toString, MetadataNumber)
      case WomBoolean(b) => new MetadataValue(b.toString, MetadataBoolean)
      case WomOptionalValue(_, Some(o)) => apply(o)
      case WomOptionalValue(_, None) => new MetadataValue("", MetadataNull)
      case value: WomValue => new MetadataValue(value.valueString, MetadataString)
      case _: Int | Long | _: java.lang.Long | _: java.lang.Integer => new MetadataValue(value.toString, MetadataInt)
      case _: Double | Float | _: java.lang.Double | _: java.lang.Float => new MetadataValue(value.toString, MetadataNumber)
      case _: Boolean | _: java.lang.Boolean => new MetadataValue(value.toString, MetadataBoolean)
      case offsetDateTime: OffsetDateTime => new MetadataValue(offsetDateTime.toUtcMilliString, MetadataString)
      case other => new MetadataValue(other.toString, MetadataString)
    }
  }
}

object MetadataType {
  val log: Logger = LoggerFactory.getLogger("Metadata Type")

  def fromString(s: String): MetadataType = s match {
    case MetadataString.typeName => MetadataString
    case MetadataInt.typeName => MetadataInt
    case MetadataNumber.typeName => MetadataNumber
    case MetadataBoolean.typeName => MetadataBoolean
    case MetadataNull.typeName => MetadataNull
    case _ =>
      log.warn(s"Unknown Metadata type $s. Falling back to MetadataString type")
      MetadataString
  }
}

final case class MetadataValue(value: String, valueType: MetadataType)

final case class MetadataEvent(key: MetadataKey, value: Option[MetadataValue], offsetDateTime: OffsetDateTime)

final case class MetadataQueryJobKey(callFqn: String, index: Option[Int], attempt: Option[Int])

object MetadataQueryJobKey {
  def forMetadataJobKey(jobKey: MetadataJobKey) = MetadataQueryJobKey(jobKey.callFqn, jobKey.index, Option(jobKey.attempt))
}

case class MetadataQuery(workflowId: WorkflowId,
                         jobKey: Option[MetadataQueryJobKey],
                         key: Option[String],
                         includeKeysOption: Option[NonEmptyList[String]],
                         excludeKeysOption: Option[NonEmptyList[String]],
                         expandSubWorkflows: Boolean)

object MetadataQuery {
  def forWorkflow(workflowId: WorkflowId) = MetadataQuery(workflowId, None, None, None, None, expandSubWorkflows = false)

  def forJob(workflowId: WorkflowId, jobKey: MetadataJobKey): MetadataQuery = {
    MetadataQuery(workflowId, Option(MetadataQueryJobKey.forMetadataJobKey(jobKey)), None, None, None, expandSubWorkflows = false)
  }

  def forKey(key: MetadataKey): MetadataQuery = {
    MetadataQuery(key.workflowId, key.jobKey map MetadataQueryJobKey.forMetadataJobKey, Option(key.key), None, None, expandSubWorkflows = false)
  }
}
