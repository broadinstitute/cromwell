package cromwell.webservice

import java.sql.Timestamp
import java.time.OffsetDateTime

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.core.{WorkflowId, WorkflowState}
import cromwell.database.obj.WorkflowMetadataKeys
import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.services.MetadataServiceActor._
import cromwell.services.ServiceRegistryActor.ServiceRegistryFailure
import cromwell.services._
import cromwell.webservice.MetadataBuilderActor.WaitingForMetadataService
import cromwell.webservice.MetadataBuilderActor.{Idle, IndexedJsonObject, IndexedJsonValue, IndexedPrimitiveJson, MetadataBuilderActorState}
import cromwell.webservice.PerRequest.{RequestComplete, RequestCompleteWithHeaders}
import org.slf4j.LoggerFactory
import spray.http.StatusCodes
import spray.httpx.SprayJsonSupport._
import spray.json.{JsString, _}
import spray.json.{JsArray, JsObject, JsValue, _}

import scala.language.postfixOps
import scala.util.Try
import scalaz.Scalaz._
import scalaz.Semigroup


object MetadataBuilderActor {
  sealed trait MetadataBuilderActorState
  case object Idle extends MetadataBuilderActorState
  case object WaitingForMetadataService extends MetadataBuilderActorState

  def props(serviceRegistryActor: ActorRef) = {
    Props(new MetadataBuilderActor(serviceRegistryActor))
  }

  val log = LoggerFactory.getLogger("MetadataBuilder")

  private val KeySeparator = ':'
  private val bracketMatcher = """\[(\d+)\]""".r
  private val startMatcher = """^([^\[]+)\[""".r
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
    }
  }

  /** Types of element supported in a dotted key notation */
  private sealed trait KeyElement {
    def toIndexedJson(value: IndexedJsonValue): IndexedJsonValue
  }

  private case class ListElement(name: String, indexes: List[String]) extends KeyElement {
    def toIndexedJson(innerValue: IndexedJsonValue) = {
      if (indexes.isEmpty) {
        IndexedJsonObject(Map(name -> innerValue))
      } else {
       /*
       * The last index is the one that the innerValue should have in the innerList.
       * From there lists are fold into one another until we reach the first index.
       * e.g l[1][2] = "a" means
       *  "l": [
       *         [
       *           "a" <- will have index 2 in the inner list
       *         ] <- inner list: will have index 1 in the outer list
       *       ]
       *
       * Important note: Indexes are used for sorting purposes ONLY
       * An index of 2 DOES NOT guarantee that a value will be at index 2 in a list
       */
        val innerList = IndexedJsonList(Vector(innerValue.withIndex(indexes.last)))
        val list = indexes.init.foldRight(innerList)((index, acc) => {
          IndexedJsonList(Vector(acc.withIndex(index)))
        })

        IndexedJsonObject(Map(name -> list))
      }
    }
  }
  private case class ObjectElement(name: String) extends KeyElement {
    def toIndexedJson(value: IndexedJsonValue) = IndexedJsonObject(Map(name -> value))
  }

  private def parseKeyChunk(chunk: String): KeyElement = {
    startMatcher.findFirstMatchIn(chunk) match {
      case Some(listNameRegex) =>
        val indexes = bracketMatcher.findAllMatchIn(chunk).map(_.group(1)).toList
        ListElement(listNameRegex.group(1), indexes)
      case _ => ObjectElement(chunk)
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
    val innerValue: IndexedJsonValue = metadataValueToIndexedJson(value)
    str.split(KeySeparator).foldRight(innerValue)((chunk, acc) => { parseKeyChunk(chunk).toIndexedJson(acc) })
  }

  private case class MetadataForAttempt(attempt: Int, metadata: IndexedJsonValue)
  /** There's one IndexedJsonValue per attempt, hence the list. */
  private case class MetadataForIndex(index: Int, metadata: List[IndexedJsonValue])

  implicit def dateTimeOrdering: Ordering[OffsetDateTime] = scala.Ordering.fromLessThan(_ isBefore _)

  implicit val timestampOrdering: Ordering[Timestamp] = scala.Ordering.fromLessThan(_.compareTo(_) < 0)

  /** Sort events by timestamp, transform them into IndexedJsonValues, and merge them together. */
  private def eventsToIndexedJson(events: Seq[MetadataEvent]): IndexedJsonValue = {
    events match {
      case empty if empty.isEmpty => IndexedJsonObject.empty
      case evts =>
        evts sortBy { _.offsetDateTime } map { e => keyValueToIndexedJson(e.key.key, e.value) } reduceLeft(_ |+| _)
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

  private def reduceWorkflowEvents(workflowEvents: Seq[MetadataEvent]): Seq[MetadataEvent] = {
    // This handles state specially so a sensible final value is returned irrespective of the order in which raw state
    // events were recorded in the journal.
    val (workflowStatusEvents, workflowNonStatusEvents) = workflowEvents partition(_.key.key == WorkflowMetadataKeys.Status)

    val ordering = implicitly[Ordering[WorkflowState]]
    // This orders by value in WorkflowState CRDT resolution, not necessarily the chronologically most recent state.
    val sortedStateEvents = workflowStatusEvents sortWith { case (a, b) => ordering.gt(a.value.toWorkflowState, b.value.toWorkflowState) }
    workflowNonStatusEvents ++ sortedStateEvents.headOption.toList
  }

  private def parseWorkflowEventsToIndexedJsonValue(events: Seq[MetadataEvent]): IndexedJsonValue = {
    // Partition if sequence of events in a pair of (Workflow level events, Call level events)
    val (workflowLevel, callLevel) = events partition { _.key.jobKey.isEmpty }
    val foldedWorkflowValues = eventsToIndexedJson(reduceWorkflowEvents(workflowLevel))

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

  /**
    * Parse a Seq of MetadataEvent into a full Json metadata response.
    */
  private def parse(events: Seq[MetadataEvent]): JsObject = {
    JsObject(events.groupBy(_.key.workflowId.toString) mapValues parseWorkflowEvents)
  }

  implicit class EnhancedMetadataValue(val value: MetadataValue) extends AnyVal {
    def toWorkflowState: WorkflowState = WorkflowState.fromString(value.value)
  }
}

class MetadataBuilderActor(serviceRegistryActor: ActorRef) extends LoggingFSM[MetadataBuilderActorState, Unit]
  with DefaultJsonProtocol with WorkflowQueryPagination {

  import WorkflowJsonSupport._

  startWith(Idle, ())
  val tag = self.path.name

  when(Idle) {
    case Event(action: MetadataServiceAction, _) =>
      serviceRegistryActor ! action
      goto(WaitingForMetadataService)
  }

  private def allDone = {
    context stop self
    stay()
  }

  when(WaitingForMetadataService) {
    case Event(MetadataLookupResponse(query, metadata), _) =>
      context.parent ! RequestComplete(StatusCodes.OK, processMetadataResponse(query, metadata))
      allDone
    case Event(StatusLookupResponse(w, status), _) =>
      context.parent ! RequestComplete(StatusCodes.OK, processStatusResponse(w, status))
      allDone
    case Event(StatusLookupNotFound(w), _) =>
      context.parent ! APIResponse.workflowNotFound(w)
      allDone
    case Event(StatusLookupFailed(_, t), _) =>
      context.parent ! RequestComplete(StatusCodes.InternalServerError, APIResponse.error(t))
      allDone
    case Event(failure: ServiceRegistryFailure, _) =>
      val response = APIResponse.fail(new RuntimeException("Can't find metadata service"))
      context.parent ! RequestComplete(StatusCodes.InternalServerError, response)
      allDone
    case Event(WorkflowQuerySuccess(uri, response, metadata), _) =>
      context.parent ! RequestCompleteWithHeaders(response, generateLinkHeaders(uri, metadata):_*)
      allDone
    case Event(WorkflowQueryFailure(t), _) =>
      context.parent ! RequestComplete(StatusCodes.InternalServerError, APIResponse.error(t))
      allDone
    case Event(unexpectedMessage, stateData) =>
      val response = APIResponse.fail(new RuntimeException(s"MetadataBuilderActor $tag(WaitingForMetadataService, $stateData) got an unexpected message: $unexpectedMessage"))
      context.parent ! RequestComplete(StatusCodes.InternalServerError, response)
      context stop self
      stay()
  }

  def processMetadataResponse(query: MetadataQuery, eventsList: Seq[MetadataEvent]): JsObject = {
    // Should we send back some message ? Or even fail the request instead ?
    if (eventsList.isEmpty) JsObject(Map.empty[String, JsValue])
    else {
      query match {
        case MetadataQuery(w, _, _) => workflowMetadataResponse(w, eventsList)
        case _ => MetadataBuilderActor.parse(eventsList)
      }
    }
  }

  def processStatusResponse(workflowId: WorkflowId, status: WorkflowState): JsObject = {
    val state: IndexedJsonValue = IndexedJsonObject(Map(WorkflowMetadataKeys.Status -> IndexedPrimitiveJson(JsString(status.toString))))
    val id: IndexedJsonValue = IndexedJsonObject(Map(WorkflowMetadataKeys.Id -> IndexedPrimitiveJson(JsString(workflowId.toString))))
    (state |+| id).toJson.asJsObject
  }

  private def workflowMetadataResponse(workflowId: WorkflowId, eventsList: Seq[MetadataEvent]) = JsObject(MetadataBuilderActor.parseWorkflowEvents(eventsList).fields + ("workflowId" -> JsString(workflowId.toString)))
}
