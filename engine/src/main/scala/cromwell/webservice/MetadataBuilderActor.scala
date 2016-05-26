package cromwell.webservice

import java.sql.Timestamp

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.core.WorkflowId
import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine.WorkflowState
import cromwell.engine.workflow.WorkflowMetadataKeys
import cromwell.services
import cromwell.services.MetadataServiceActor._
import cromwell.services.MetadataValue._
import cromwell.services.ServiceRegistryActor.ServiceRegistryFailure
import cromwell.services._
import cromwell.webservice.MetadataBuilderActor.{Idle, IndexedJsonObject, MetadataBuilderActorState, WaitingForMetadata}
import cromwell.webservice.PerRequest.RequestComplete
import org.joda.time.DateTime
import org.slf4j.LoggerFactory
import spray.http.StatusCodes
import spray.httpx.SprayJsonSupport._
import spray.json._

import scala.annotation.tailrec
import scala.language.postfixOps
import scala.util.{Failure, Try}
import scalaz.Scalaz._
import scalaz.Semigroup

object MetadataBuilderActor {
  sealed trait MetadataBuilderActorState
  case object Idle extends MetadataBuilderActorState
  case object WaitingForMetadata extends MetadataBuilderActorState

  def props(serviceRegistryActor: ActorRef) = {
    Props(new MetadataBuilderActor(serviceRegistryActor))
  }

  val log = LoggerFactory.getLogger("MetadataBuilder")

  private val KeySeparator = ':'
  private val ListNameCaptureID = "listName"
  private val ListIndexCaptureID = "listIndex"
  private val ListMatcher = "^(.+)\\[([a-zA-Z0-9]+)\\]$".r(ListNameCaptureID, ListIndexCaptureID)
  private val AttemptKey = "attempt"
  private val ShardKey = "shardIndex"

  // Define how 2 IndexedJsonValues can be summed
  private implicit val IndexedJsonSemigroup: Semigroup[IndexedJsonValue] = {
    new Semigroup[IndexedJsonValue] {
      def append(f1: IndexedJsonValue, f2: => IndexedJsonValue): IndexedJsonValue = {
        (f1, f2) match {
          case (o1: IndexedJsonObject, o2: IndexedJsonObject) =>
            val sg = implicitly[Semigroup[Map[String, IndexedJsonValue]]]
            IndexedJsonObject(sg.append(o1.v, o2.v), o1.index)
          case (o1: IndexedJsonList, o2: IndexedJsonList) =>
            val sg = implicitly[Semigroup[IndexedJsonValue]]
            // Try to merge together values with the same index
            val merged = (o1.v ++ o2.v).groupBy(_.index) map { case (_, v) => v.reduceLeft(sg.append(_, _)) }
            val sorted = merged.toVector sortBy { _.index }
            IndexedJsonList(sorted, o1.index)
          case (o1, o2) => o2 // We assume values are sorted by timestamp, so if we can't merge we keep the new one.
        }
      }
    }
  }

  /** Types of element supported in a dotted key notation */
  private sealed trait KeyElement
  private case class ListElement(name: String, index: String) extends KeyElement
  private case class ObjectElement(name: String) extends KeyElement

  private def parseKeyChunk(chunk: String): KeyElement = {
    ListMatcher.findAllIn(chunk) match {
      case matchIterator if matchIterator.hasNext => ListElement(matchIterator.group(ListNameCaptureID), matchIterator.group(ListIndexCaptureID))
      case _ => ObjectElement(chunk)
    }
  }

  /** Simplified version of Json data structure
    * Every value has an index property that enables late sorting of arrays. */
  private sealed trait IndexedJsonValue {
    def index:String
    def toJson: JsValue
    def withIndex(index: String): IndexedJsonValue
  }
  private case class IndexedJsonList(v: Vector[IndexedJsonValue], index: String = "") extends IndexedJsonValue {
    override val toJson = JsArray(v map { _.toJson })
    override def withIndex(index: String) = this.copy(index = index)
  }
  private object IndexedJsonObject {
    def empty = IndexedJsonObject(Map.empty[String, IndexedJsonValue])
  }
  private case class IndexedJsonObject(v: Map[String, IndexedJsonValue], index: String = "") extends IndexedJsonValue {
    override val toJson = JsObject(v mapValues { _.toJson })
    override def withIndex(index: String) = this.copy(index = index)
  }
  private case class IndexedPrimitiveJson(v: JsValue, index: String = "") extends IndexedJsonValue {
    override val toJson = v
    override def withIndex(index: String) = this.copy(index = index)
  }

  /**
    * Create an IndexedJson structure from a list of key chunks.
    */
  @tailrec
  private def parseKeyRec(parts: List[String], current: IndexedJsonValue): IndexedJsonValue = {
    parts match {
      case p :: r =>
        parseKeyChunk(p) match {
          case ObjectElement(name) => parseKeyRec(r, IndexedJsonObject(Map(name -> current)))
          case ListElement(name, index) => parseKeyRec(r, IndexedJsonObject(Map(name -> IndexedJsonList(Vector(current.withIndex(index))))))
        }
      case Nil => current
    }
  }

  private def metadataValueToIndexedJson(value: MetadataValue): IndexedPrimitiveJson = {
    val coerced = value.valueType match {
      case MetadataInt => Try(IndexedPrimitiveJson(JsNumber(value.value.toInt)))
      case MetadataNumber => Try(IndexedPrimitiveJson(JsNumber(value.value.toDouble)))
      case MetadataBoolean => Try(IndexedPrimitiveJson(JsBoolean(value.value.toBoolean)))
      case MetadataString => Try(IndexedPrimitiveJson(JsString(value.value)))
    }

    coerced recover {
      case e =>
      log.warn(s"Failed to coerce ${value.value} to ${value.valueType}. Falling back to String.", e)
      IndexedPrimitiveJson(JsString(value.value))
    } get
  }

  private def keyValueToIndexedJson(str: String, value: MetadataValue): IndexedJsonValue = {
    // The list is reversed because we build the IndexedJson from the bottom up
    val reversed = str.split(KeySeparator).toList.reverse
    parseKeyRec(reversed, metadataValueToIndexedJson(value))
  }

  private case class MetadataForAttempt(attempt: Int, metadata: IndexedJsonValue)
  /** There's one IndexedJsonValue per attempt, hence the list. */
  private case class MetadataForIndex(index: Int, metadata: List[IndexedJsonValue])

  implicit def dateTimeOrdering: Ordering[DateTime] = scala.Ordering.fromLessThan(_ isBefore _)

  implicit val timestampOrdering: Ordering[Timestamp] = scala.Ordering.fromLessThan(_.compareTo(_) < 0)

  /** Sort events by timestamp, transform them into IndexedJsonValues, and merge them together. */
  private def eventsToIndexedJson(events: Seq[MetadataEvent]): IndexedJsonValue = {
    events match {
      case empty if empty.isEmpty => IndexedJsonObject.empty
      case evts => evts sortBy { _.timestamp } map { e => keyValueToIndexedJson(e.key.key, e.value) } reduceLeft(_ |+| _)
    }
  }

  private def eventsToAttemptMetadata(attempt: Int, events: Seq[MetadataEvent]) = {
    val withAttemptField = eventsToIndexedJson(events) |+| IndexedJsonObject(Map(AttemptKey -> IndexedPrimitiveJson(JsNumber(attempt))))
    MetadataForAttempt(attempt, withAttemptField)
  }

  private def attemptMetadataToIndexMetadata(index: ExecutionIndex, attemptMetadata: Iterable[MetadataForAttempt]) = {
    def addIndexProperty(value: IndexedJsonValue) = value |+| IndexedJsonObject(Map(ShardKey -> IndexedPrimitiveJson(JsNumber(index.getOrElse(-1)))))
    val metadata = attemptMetadata.toList.sortBy(_.attempt) map { mdForAttempt => addIndexProperty(mdForAttempt.metadata) }
    MetadataForIndex(index.getOrElse(-1), metadata)
  }

  private def parseWorkflowEventsToIndexedJsonValue(events: Seq[MetadataEvent]): IndexedJsonValue = {
    // Partition if sequence of events in a pair of (Workflow level events, Call level events)
    val (workflowLevel, callLevel) = events partition { _.key.jobKey.isEmpty }
    val foldedWorkflowValues = eventsToIndexedJson(workflowLevel)

    val callsGroupedByFQN = callLevel groupBy { _.key.jobKey.get.callFqn }
    val callsGroupedByFQNAndIndex = callsGroupedByFQN mapValues { _ groupBy { _.key.jobKey.get.index } }
    val callsGroupedByFQNAndIndexAndAttempt = callsGroupedByFQNAndIndex mapValues { _ mapValues { _ groupBy { _.key.jobKey.get.attempt } } }

    val callsMap = callsGroupedByFQNAndIndexAndAttempt mapValues { eventsForIndex =>
      eventsForIndex mapValues { eventsForAttempt =>
        eventsForAttempt map Function.tupled(eventsToAttemptMetadata)
      } map { Function.tupled(attemptMetadataToIndexMetadata) }
    } mapValues { md => IndexedJsonList(md.toVector.sortBy(_.index) flatMap { _.metadata }) }

    val callsObject = IndexedJsonObject(callsMap)
    val wrappedCalls = IndexedJsonObject(Map(WorkflowMetadataKeys.Calls -> callsObject))

    foldedWorkflowValues |+| wrappedCalls
  }

  private def parseWorkflowEvents(events: Seq[MetadataEvent]): JsObject = parseWorkflowEventsToIndexedJsonValue(events).toJson.asJsObject

  private def foldStates(eventsList: Seq[MetadataEvent]): IndexedJsonValue = {
    import WorkflowState._
    val state: WorkflowState = eventsList collect { case MetadataEvent(_, v, _) => WorkflowState.fromString(v.value) } reduceLeft(_ |+| _)
    IndexedJsonObject(Map(WorkflowMetadataKeys.Status -> IndexedPrimitiveJson(JsString(state.toString))))
  }

  private def foldStatesAndAddWorkflowId(workflowId: WorkflowId, eventsList: Seq[MetadataEvent]): JsObject = {
    val result = foldStates(eventsList) |+| IndexedJsonObject(Map(WorkflowMetadataKeys.Id -> IndexedPrimitiveJson(JsString(workflowId.toString))))
    result.toJson.asJsObject
  }

  /**
    * Parse a Seq of MetadataEvent into a full Json metadata response.
    */
  private def parse(events: Seq[MetadataEvent]): JsObject = {
    JsObject(events.groupBy(_.key.workflowId.toString) mapValues parseWorkflowEvents)
  }

}

class MetadataBuilderActor(serviceRegistryActor: ActorRef) extends LoggingFSM[MetadataBuilderActorState, Unit] with DefaultJsonProtocol {
  // Don't believe IntelliJ, this is used.
  import WorkflowJsonSupport._

  startWith(Idle, Unit)

  when(Idle) {
    case Event(action: MetadataServiceAction, _) =>
      serviceRegistryActor ! action
      goto(WaitingForMetadata)
  }

  when(WaitingForMetadata) {
    case Event(MetadataLookupResponse(query, metadata), _) =>
      context.parent ! RequestComplete(StatusCodes.OK, processMetadataResponse(query, metadata))
      context stop self
      stay()
    case Event(failure: ServiceRegistryFailure, _) =>
      val response = APIResponse.fail(new RuntimeException("Can't find metadata service"))
      context.parent ! RequestComplete(StatusCodes.InternalServerError, response)
      context stop self
      stay()
  }

  def processMetadataResponse(query: MetadataQuery, eventsList: Seq[MetadataEvent]): JsObject = {
    // Should we send back some message ? Or even fail the request instead ?
    if (eventsList.isEmpty) JsObject(Map.empty[String, JsValue])
    else {
      query match {
        case MetadataQuery(w, None, None) => IndexedJsonObject(Map(w.id.toString -> MetadataBuilderActor.parseWorkflowEventsToIndexedJsonValue(eventsList))).toJson.asJsObject
        case MetadataQuery(w, None, Some(WorkflowMetadataKeys.Status)) => MetadataBuilderActor.foldStatesAndAddWorkflowId(w, eventsList)
        case _ => MetadataBuilderActor.parse(eventsList)
      }
    }
  }
}
