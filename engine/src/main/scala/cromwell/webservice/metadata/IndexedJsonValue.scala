package cromwell.webservice.metadata

import spray.json._

import scalaz.{Monoid, Semigroup}
import scalaz.Scalaz._

private object IndexedJsonValue {
   implicit val IndexedJsonSemigroup: Semigroup[IndexedJsonValue] = {
    new Monoid[IndexedJsonValue] {
      def append(f1: IndexedJsonValue, f2: => IndexedJsonValue): IndexedJsonValue = {
        (f1, f2) match {
          case (o1: IndexedJsonObject, o2: IndexedJsonObject) =>
            val sg = implicitly[Semigroup[Map[String, IndexedJsonValue]]]
            IndexedJsonObject(sg.append(o1.v, o2.v), o1.index)
          case (o1: IndexedJsonList, o2: IndexedJsonList) =>
            val sg = implicitly[Semigroup[IndexedJsonValue]]
            /*
             * Try to merge together values with the same index (if defined)
             *
             * If no index is defined the values are just left as is
             * Note that as long as we require the index to be non empty for a list element (ie disallow l[]), this should not happen
             */
            val merged = (o1.v ++ o2.v).groupBy(_.index) flatMap {
              case (index, v) if index.isDefined => Vector(v.reduceLeft(sg.append(_, _)))
              case (_, v) => v
            }
            val sorted = merged.toVector sortBy { _.index }
            IndexedJsonList(sorted, o1.index)
          case (o1, o2) => o2 // We assume values are sorted by timestamp, so if we can't merge we keep the new one.
        }
      }

      override def zero: IndexedJsonValue = IndexedJsonObject.empty
    }
  }
}

/** Simplified version of Json data structure
  * Every value has an index property that enables late sorting of arrays. */
private sealed trait IndexedJsonValue {
  def index: Option[Int]
  def toJson: JsValue
  def withIndex(index: String): IndexedJsonValue

  // The regex guarantees that indexes contain only digits, so toInt will succeed as long as the index is < Int.max
  def toIndex(str: String) = str match {
    case s if s.isEmpty => None
    case s => Option(s.toInt)
  }
}

private case class IndexedEmptyJson(index: Option[Int] = None) extends IndexedJsonValue {
  override def withIndex(index: String): IndexedJsonValue = this.copy(index = toIndex(index))
  override def toJson: JsValue = JsObject(Map.empty[String, JsValue])
}

private object IndexedJsonList {
  def empty = IndexedJsonList(Vector.empty[IndexedJsonValue])
}

private case class IndexedJsonList(v: Vector[IndexedJsonValue], index: Option[Int] = None) extends IndexedJsonValue {
  override val toJson = JsArray(v map { _.toJson })
  override def withIndex(index: String) = this.copy(index = toIndex(index))
}

private object IndexedJsonObject {
  def empty = IndexedJsonObject(Map.empty[String, IndexedJsonValue])
}

private case class IndexedJsonObject(v: Map[String, IndexedJsonValue], index: Option[Int] = None) extends IndexedJsonValue {
  override val toJson = JsObject(v mapValues { _.toJson })
  override def withIndex(index: String) = this.copy(index = toIndex(index))
}

private case class IndexedPrimitiveJson(v: JsValue, index: Option[Int] = None) extends IndexedJsonValue {
  override val toJson = v
  override def withIndex(index: String) = this.copy(index = toIndex(index))
}