package cromwell.webservice.metadata

import java.time.OffsetDateTime

import cats.{Monoid, Semigroup}
import cats.instances.map._
import cromwell.services.metadata.CallMetadataKeys
import spray.json._


object IndexedJsonValue {
  private implicit val dateTimeOrdering: Ordering[OffsetDateTime] = scala.Ordering.fromLessThan(_ isBefore _)
  private val timestampedJsValueOrdering: Ordering[TimestampedJsValue] = scala.Ordering.by(_.timestamp)

  implicit val TimestampedJsonMonoid: Monoid[TimestampedJsValue] = new Monoid[TimestampedJsValue] {
    def combine(f1: TimestampedJsValue, f2: TimestampedJsValue): TimestampedJsValue = {
      (f1, f2) match {
        case (o1: TimestampedJsObject, o2: TimestampedJsObject) =>
          val sg = implicitly[Semigroup[Map[String, TimestampedJsValue]]]
          TimestampedJsObject(sg.combine(o1.v, o2.v), dateTimeOrdering.max(o1.timestamp, o2.timestamp))
        case (o1: TimestampedJsList, o2: TimestampedJsList) =>
          val sg = implicitly[Semigroup[Map[Int, TimestampedJsValue]]]
          TimestampedJsList(sg.combine(o1.v, o2.v), dateTimeOrdering.max(o1.timestamp, o2.timestamp))
        case (o1, o2) => timestampedJsValueOrdering.max(o1, o2)
      }
    }

    override def empty: TimestampedJsValue = TimestampedJsObject(Map.empty, OffsetDateTime.now)
  }
}

/** Customized version of Json data structure, to account for timestamped values and lazy array creation */
sealed trait TimestampedJsValue {
  def toJson(expandedValues: Map[String, JsValue]): JsValue
  def timestamp: OffsetDateTime
}

private case class TimestampedJsList(v: Map[Int, TimestampedJsValue], timestamp: OffsetDateTime) extends TimestampedJsValue {
  override def toJson(expandedValues: Map[String, JsValue]) = JsArray(v.values.toVector map { _.toJson(expandedValues) })
}

private case class TimestampedJsObject(v: Map[String, TimestampedJsValue], timestamp: OffsetDateTime) extends TimestampedJsValue {
  override def toJson(expandedValues: Map[String, JsValue]) = {
    val mappedValues = v map {
      case (key, subWorkflowId: TimestampedJsPrimitive) if key == CallMetadataKeys.SubWorkflowId =>
        val subId = subWorkflowId.v.asInstanceOf[JsString]
        expandedValues.get(subId.value) map { subMetadata =>
           CallMetadataKeys.SubWorkflowMetadata -> subMetadata
        } getOrElse {
          key -> subWorkflowId.v
        }
      case (key, value) => key -> value.toJson(expandedValues)
    }
    
    JsObject(mappedValues)
  }
}

private class TimestampedJsPrimitive(val v: JsValue, val timestamp: OffsetDateTime) extends TimestampedJsValue {
  override def toJson(expandedValues: Map[String, JsValue]) = v
}

private case class TimestampedEmptyJson(override val timestamp: OffsetDateTime) extends TimestampedJsPrimitive(JsObject(Map.empty[String, JsValue]), timestamp)