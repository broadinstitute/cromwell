package cromwell.webservice.metadata

import cats.{Monoid, Semigroup}
import cats.instances.map._
import cromwell.core.{ExecutionStatus, WorkflowState}
import spray.json.{JsArray, _}

object MetadataComponent {
  implicit val MetadataComponentMonoid: Monoid[MetadataComponent] = new Monoid[MetadataComponent] {
    private lazy val stringValueMapSg = implicitly[Semigroup[Map[String, MetadataComponent]]]
    private lazy val intValueMapSg = implicitly[Semigroup[Map[Int, MetadataComponent]]]
    
    def combine(f1: MetadataComponent, f2: MetadataComponent): MetadataComponent = {
      (f1, f2) match {
        case (MetadataObject(v1), MetadataObject(v2)) => MetadataObject(stringValueMapSg.combine(v1, v2))
        case (MetadataList(v1), MetadataList(v2)) => MetadataList(intValueMapSg.combine(v1, v2))
          // If there's a custom ordering, use it
        case (v1 @ MetadataPrimitive(_, Some(o1)), v2 @ MetadataPrimitive(_, Some(o2))) if o1 == o2 => o1.max(v1, v2)
          // Otherwise assume it's ordered by default and take the new one
        case (o1, o2) => o2
      }
    }

    override def empty: MetadataComponent = MetadataObject.empty
  }
  
  implicit val jsonWriter: JsonWriter[MetadataComponent] = JsonWriter.func2Writer[MetadataComponent] {
    case MetadataList(values) => JsArray(values.values.toVector map { _.toJson(this.jsonWriter) })
    case MetadataObject(values) => JsObject(values.mapValues(_.toJson(this.jsonWriter)))
    case MetadataPrimitive(value, _) => value
    case MetadataEmpty => JsObject.empty
  }
}

sealed trait MetadataComponent
case object MetadataEmpty extends MetadataComponent
object MetadataObject { def empty = MetadataObject(Map.empty) }
case class MetadataObject(v: Map[String, MetadataComponent]) extends MetadataComponent
case class MetadataList(v: Map[Int, MetadataComponent]) extends MetadataComponent

object MetadataPrimitive {
  val ExecutionStatusOrdering: Ordering[MetadataPrimitive] = Ordering.by { primitive: MetadataPrimitive =>
    ExecutionStatus.withName(primitive.v.asInstanceOf[JsString].value)
  }

  val WorkflowStateOrdering: Ordering[MetadataPrimitive] = Ordering.by { primitive: MetadataPrimitive =>
    WorkflowState.fromString(primitive.v.asInstanceOf[JsString].value)
  }
}
case class MetadataPrimitive(v: JsValue, customOrdering: Option[Ordering[MetadataPrimitive]] = None) extends MetadataComponent
