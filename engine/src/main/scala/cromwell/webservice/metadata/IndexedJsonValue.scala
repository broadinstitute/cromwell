package cromwell.webservice.metadata

import java.time.OffsetDateTime

import cats.{Monoid, Semigroup}
import cats.instances.map._
import spray.json._


private object IndexedJsonValue {
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
private sealed trait TimestampedJsValue {
  def toJson: JsValue
  def timestamp: OffsetDateTime
}

private case class TimestampedJsList(v: Map[Int, TimestampedJsValue], timestamp: OffsetDateTime) extends TimestampedJsValue {
  override val toJson = JsArray(v.values.toVector map { _.toJson })
}

private case class TimestampedJsObject(v: Map[String, TimestampedJsValue], timestamp: OffsetDateTime) extends TimestampedJsValue {
  override val toJson = JsObject(v mapValues { _.toJson })
}

private class TimestampedJsPrimitive(val v: JsValue, val timestamp: OffsetDateTime) extends TimestampedJsValue {
  override val toJson = v
}

private case class TimestampedEmptyJson(override val timestamp: OffsetDateTime) extends TimestampedJsPrimitive(JsObject(Map.empty[String, JsValue]), timestamp)