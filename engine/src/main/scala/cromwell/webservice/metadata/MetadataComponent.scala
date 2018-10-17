package cromwell.webservice.metadata

import cats.instances.map._
import cats.{Monoid, Semigroup}
import common.collections.EnhancedCollections._
import cromwell.core._
import cromwell.core.simpleton.WomValueSimpleton._
import cromwell.services.metadata._
import spray.json._

import scala.collection.immutable.TreeMap
import scala.concurrent.ExecutionContext
import scala.util.{Random, Try}

object MetadataComponent {
  private val CallsKey = "calls"
  private val AttemptKey = "attempt"
  private val ShardKey = "shardIndex"
  
  implicit val MetadataComponentMonoid: Monoid[MetadataComponent] = new Monoid[MetadataComponent] {
    private lazy val stringKeyMapSg = implicitly[Semigroup[Map[String, MetadataComponent]]]
    private lazy val intKeyMapSg = implicitly[Semigroup[Map[Int, MetadataComponent]]]

    def combine(f1: MetadataComponent, f2: MetadataComponent): MetadataComponent = {
      (f1, f2) match {
        case (MetadataObject(v1), MetadataObject(v2)) => MetadataObject(stringKeyMapSg.combine(v1, v2))
        case (MetadataList(v1), MetadataList(v2)) => MetadataList(intKeyMapSg.combine(v1, v2))
        // If there's a custom ordering, use it
        case (v1 @ MetadataPrimitive(_, Some(o1)), v2 @ MetadataPrimitive(_, Some(o2))) if o1 == o2 => o1.max(v1, v2)
        // Otherwise assume it's ordered by default and take the new one
        case (_, o2) => o2
      }
    }

    override def empty: MetadataComponent = MetadataObject.empty
  }

  val metadataPrimitiveJsonWriter: JsonWriter[MetadataPrimitive] = JsonWriter.func2Writer[MetadataPrimitive] {
    case MetadataPrimitive(MetadataValue(value, MetadataInt), _) => Try(value.toInt) map JsNumber.apply getOrElse JsString(value)
    case MetadataPrimitive(MetadataValue(value, MetadataNumber), _) => Try(value.toDouble) map JsNumber.apply getOrElse JsString(value)
    case MetadataPrimitive(MetadataValue(value, MetadataBoolean), _) => Try(value.toBoolean) map JsBoolean.apply getOrElse JsString(value)
    case MetadataPrimitive(MetadataValue(value, MetadataString), _) => JsString(value)
    case MetadataPrimitive(MetadataValue(_, MetadataNull), _) => JsNull
  }

  /* ******************************* */
  /* *** Metadata Events Parsing *** */
  /* ******************************* */

  private val KeySeparator = MetadataKey.KeySeparator
  // Split on every unescaped KeySeparator
  val KeySplitter = s"(?<!\\\\)$KeySeparator"
  private val bracketMatcher = """\[(\d*)\]""".r

  private def parseKeyChunk(chunk: String, innerValue: MetadataComponent): MetadataObject = {
    chunk.indexOf('[') match {
      // If there's no bracket, it's an object. e.g.: "calls"
      case -1 => MetadataObject(Map(chunk -> innerValue))
      // If there's a bracket it's a named list. e.g.: "executionEvents[0][1]"
      case bracketIndex =>
        // Name: "executionEvents"
        val objectName = chunk.substring(0, bracketIndex)

        // Brackets: "[0][1]"
        val brackets = chunk.substring(bracketIndex)
        // Indices as a list: List(0, 1)
        val listIndices = for {
          m <- bracketMatcher.findAllMatchIn(brackets)
          // It's possible for a bracket pair to be empty, in which case we just give it a random number
          asInt = if (m.group(1).isEmpty) Random.nextInt() else m.group(1).toInt
        } yield asInt
        // Fold into a MetadataList: MetadataList(0 -> MetadataList(1 -> innerValue))
        val init = if (innerValue == MetadataEmptyComponent) {
          // Empty value means empty list
          listIndices.toList.init.foldRight(MetadataList.empty)((index, acc) => MetadataList(TreeMap(index -> acc)))
        } else {
          listIndices.toList.foldRight(innerValue)((index, acc) => MetadataList(TreeMap(index -> acc)))
        }
        val metadataList = listIndices.toList.foldRight(init)((index, acc) => MetadataList(TreeMap(index -> acc)))

        MetadataObject(Map(objectName -> metadataList))
    }
  }

  private def customOrdering(event: MetadataEvent): Option[Ordering[MetadataPrimitive]] = event match {
    case MetadataEvent(MetadataKey(_, Some(_), key), _, _) if key == CallMetadataKeys.ExecutionStatus => Option(MetadataPrimitive.ExecutionStatusOrdering)
    case MetadataEvent(MetadataKey(_, None, key), _, _) if key == WorkflowMetadataKeys.Status => Option(MetadataPrimitive.WorkflowStateOrdering)
    case _ => None
  }

  def toMetadataComponent(query: MetadataQuery, ec: ExecutionContext)(event: MetadataEvent) = {
    lazy val primitive = event.value map { MetadataPrimitive(_, customOrdering(event)) } getOrElse MetadataEmptyComponent

    // If the event is a sub workflow id event, we might need to replace the id with a placeholder for the sub workflow metadata
    lazy val subWorkflowKey = event.key.key.replace(CallMetadataKeys.SubWorkflowId, CallMetadataKeys.SubWorkflowMetadata)
    lazy val subWorkflowComponent = event.extractSubWorkflowId.map(MetadataPlaceholderComponent.apply).getOrElse(primitive)

    val (key, component): (String, MetadataComponent) = if (event.isSubWorkflowId) 
      subWorkflowKey -> subWorkflowComponent 
    else
      event.key.key -> primitive

    contextualize(event.key, fromMetadataKeyAndPrimitive(key, component).asInstanceOf[MetadataObject])
  }

  def toMetadataComponent(event: MetadataEvent): MetadataComponent = {
    lazy val primitive = event.value map { MetadataPrimitive(_, customOrdering(event)) } getOrElse MetadataEmptyComponent
    contextualize(event.key, fromMetadataKeyAndPrimitive(event.key.key, primitive).asInstanceOf[MetadataObject])
  }

  /** Sort events by timestamp, transform them into MetadataComponent, and merge them together. */
  def apply(events: Seq[MetadataEvent], subWorkflowMetadata: Map[String, JsValue] = Map.empty): MetadataComponent = {
    MetadataEmptyComponent
  }

  def fromMetadataKeyAndPrimitive(metadataKey: String, innerComponent: MetadataComponent): MetadataComponent = {
    metadataKey.split(KeySplitter).map(_.unescapeMeta).foldRight(innerComponent)(parseKeyChunk)
  }

  private def contextualize(key: MetadataKey, component: MetadataObject): MetadataComponent = {
    key.jobKey match {
      case None => component
      case Some(jobKey) =>
        MetadataObject(
          CallsKey -> MetadataObject(
            jobKey.callFqn -> MetadataList(
              TreeMap(
                jobKey.index.getOrElse(-1) ->
                  component.copy(v = component.v ++ Map(
                    ShardKey -> MetadataPrimitive(MetadataValue(jobKey.index.getOrElse(-1))),
                    AttemptKey -> MetadataPrimitive(MetadataValue(jobKey.attempt))
                  ))
              )
            )
          )
        )
    }
  }
  
  implicit class EnhancedMetadataEvent(val event: MetadataEvent) extends AnyVal {
    def isSubWorkflowId: Boolean = event.key.key.endsWith(CallMetadataKeys.SubWorkflowId)
    def extractSubWorkflowId: Option[WorkflowId] = {
      event.value.map(value => WorkflowId.fromString(value.value))
    }
  }
  
  implicit val DefaultMetadataComponentJsonWriter = new MetadataComponentJsonWriter(Map.empty)
  
  class MetadataComponentJsonWriter(subWorkflowsMapping: Map[WorkflowId, JsObject]) extends JsonWriter[MetadataComponent] {
    override def write(obj: MetadataComponent) = obj match {
      case MetadataList(values) => JsArray(values.values.toVector map { _.toJson(this) })
      case MetadataObject(values) => JsObject(values.safeMapValues(_.toJson(this)))
      case primitive: MetadataPrimitive => metadataPrimitiveJsonWriter.write(primitive)
      case MetadataEmptyComponent => JsObject.empty
      case MetadataNullComponent => JsNull
      case MetadataPlaceholderComponent(subWorkflowId) => subWorkflowsMapping.getOrElse(subWorkflowId, JsNull)
      case MetadataJsonComponent(jsValue) => jsValue
    }
  }
}

sealed trait MetadataComponent
case object MetadataEmptyComponent extends MetadataComponent
case object MetadataNullComponent extends MetadataComponent
case class MetadataPlaceholderComponent(subWorkflowId: WorkflowId) extends MetadataComponent

// Metadata Object  
object MetadataObject {
  def empty = new MetadataObject(Map.empty)
  def apply(kvPair: (String, MetadataComponent)*) = {
    new MetadataObject(kvPair.toMap)
  }
}

case class MetadataObject(v: Map[String, MetadataComponent]) extends MetadataComponent

// Metadata List
object MetadataList {
  def empty = new MetadataList(Map.empty)
  def apply(components: List[MetadataComponent]) = new MetadataList(components.zipWithIndex.map({case (c, i) => i -> c}).toMap)
}
case class MetadataList(v: Map[Int, MetadataComponent]) extends MetadataComponent

// Metadata Primitive
object MetadataPrimitive {
  val ExecutionStatusOrdering: Ordering[MetadataPrimitive] = Ordering.by { primitive: MetadataPrimitive =>
    ExecutionStatus.withName(primitive.v.value)
  }

  val WorkflowStateOrdering: Ordering[MetadataPrimitive] = Ordering.by { primitive: MetadataPrimitive =>
    WorkflowState.withName(primitive.v.value)
  }
}
case class MetadataPrimitive(v: MetadataValue, customOrdering: Option[Ordering[MetadataPrimitive]] = None) extends MetadataComponent

// Metadata Component that owns an already computed JsValue
case class MetadataJsonComponent(jsValue: JsValue) extends MetadataComponent
