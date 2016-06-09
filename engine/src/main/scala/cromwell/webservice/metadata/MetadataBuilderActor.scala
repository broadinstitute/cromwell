package cromwell.webservice.metadata

import java.time.OffsetDateTime

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.core.{ExecutionIndex, WorkflowId, WorkflowState}
import cromwell.database.obj.WorkflowMetadataKeys
import ExecutionIndex.ExecutionIndex
import cromwell.services.MetadataServiceActor._
import cromwell.services.ServiceRegistryActor.ServiceRegistryFailure
import cromwell.services._
import cromwell.webservice.PerRequest.{RequestComplete, RequestCompleteWithHeaders}
import cromwell.webservice.metadata.MetadataBuilderActor.{Idle, WaitingForMetadataService, MetadataBuilderActorState}
import cromwell.webservice.{WorkflowJsonSupport, APIResponse, WorkflowQueryPagination}
import org.slf4j.LoggerFactory
import spray.http.StatusCodes
import spray.httpx.SprayJsonSupport._
import spray.json._
import IndexedJsonValue._

import scala.collection.immutable.TreeMap
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import scalaz.Scalaz._

object MetadataBuilderActor {
  sealed trait MetadataBuilderActorState
  case object Idle extends MetadataBuilderActorState
  case object WaitingForMetadataService extends MetadataBuilderActorState

  def props(serviceRegistryActor: ActorRef) = {
    Props(new MetadataBuilderActor(serviceRegistryActor))
  }

  val log = LoggerFactory.getLogger("MetadataBuilder")

  private val KeySeparator = ':'
  private val bracketMatcher = """\[(\d*)\]""".r
  private val startMatcher = """^([^\[]+)\[""".r
  private val AttemptKey = "attempt"
  private val ShardKey = "shardIndex"

  /** Types of element supported in a dotted key notation */
  private sealed trait KeyElement {
    def toIndexedJson(value: TimestampedJsValue): TimestampedJsValue
  }

  private case class ListElement(name: String, indexes: List[String]) extends KeyElement {
    def toIndexedJson(innerValue: TimestampedJsValue) = {
      if (indexes.isEmpty) {
        TimestampedJsObject(Map(name -> innerValue), innerValue.timestamp)
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
        val list = innerValue match {
           // Empty value in a list means empty list
          case TimestampedEmptyJson(timestamp) => TimestampedJsList(Map.empty, timestamp)
          case nonEmptyValue =>
             /*
              * This creates a (possibly nested) list, by folding over the indexes.
              * The resulting list will be as deep as there are elements in "indexes"
              * First we create the deepest list, that will contain innerValue (the actual value we want in the list)
              * e.g with l[1][2] = "a". indexes will be List[1, 2]. innerValue will be "a".
              * innerList is TimestampedJsList(Map(2 -> "a"), [timestamp of a])
              */
            val innerList = TimestampedJsList(TreeMap(indexes.last.toInt -> innerValue), nonEmptyValue.timestamp)
            /* Then, stating with this innerList, we wrap around it as many lists as (indexes.length - 1) (because we used the last index for the innerValue above)
             * Continuing with this example, result will be TimestampedJsList(Map(1 -> TimestampedJsList(Map(2 -> "a"))))
             */
            indexes.init.foldRight(innerList)((index, acc) => {
              TimestampedJsList(TreeMap(index.toInt -> acc), acc.timestamp)
            })
        }

        TimestampedJsObject(Map(name -> list), list.timestamp)
      }
    }
  }
  private case class ObjectElement(name: String) extends KeyElement {
    def toIndexedJson(value: TimestampedJsValue) = TimestampedJsObject(Map(name -> value), value.timestamp)
  }

  private def parseKeyChunk(chunk: String): KeyElement = {
    startMatcher.findFirstMatchIn(chunk) match {
      case Some(listNameRegex) =>
        val indexes = bracketMatcher.findAllMatchIn(chunk).map(_.group(1)).toList
        ListElement(listNameRegex.group(1), indexes)
      case _ => ObjectElement(chunk)
    }
  }

  private def metadataValueToIndexedJson(value: Option[MetadataValue], timestamp: OffsetDateTime): TimestampedJsValue = {
    value map { someValue =>
      val coerced: Try[TimestampedJsPrimitive] = someValue.valueType match {
        case MetadataInt => Try(new TimestampedJsPrimitive(JsNumber(someValue.value.toInt), timestamp))
        case MetadataNumber => Try(new TimestampedJsPrimitive(JsNumber(someValue.value.toDouble), timestamp))
        case MetadataBoolean => Try(new TimestampedJsPrimitive(JsBoolean(someValue.value.toBoolean), timestamp))
        case MetadataString => Try(new TimestampedJsPrimitive(JsString(someValue.value), timestamp))
      }

      coerced match {
        case Success(v) => v
        case Failure(e) =>
          log.warn(s"Failed to coerce ${someValue.value} to ${someValue.valueType}. Falling back to String.", e)
          new TimestampedJsPrimitive(JsString(someValue.value), timestamp)
      }
    } getOrElse TimestampedEmptyJson(timestamp)
  }

  private def keyValueToIndexedJson(str: String, value: Option[MetadataValue], timestamp: OffsetDateTime): TimestampedJsValue = {
    val innerValue: TimestampedJsValue = metadataValueToIndexedJson(value, timestamp)
    str.split(KeySeparator).foldRight(innerValue)((chunk, acc) => { parseKeyChunk(chunk).toIndexedJson(acc) })
  }

  private case class MetadataForAttempt(attempt: Int, metadata: JsObject)
  /** There's one TimestampedJsValue per attempt, hence the list. */
  private case class MetadataForIndex(index: Int, metadata: List[JsObject])

  implicit val dateTimeOrdering: Ordering[OffsetDateTime] = scala.Ordering.fromLessThan(_ isBefore _)

  /** Sort events by timestamp, transform them into TimestampedJsValues, and merge them together. */
  private def eventsToIndexedJson(events: Seq[MetadataEvent]): TimestampedJsValue = {
     events map { e => keyValueToIndexedJson(e.key.key, e.value, e.offsetDateTime) } reduce(_ |+| _)
  }

  private def eventsToAttemptMetadata(attempt: Int, events: Seq[MetadataEvent]) = {
    val withAttemptField = JsObject(eventsToIndexedJson(events).toJson.asJsObject.fields + (AttemptKey -> JsNumber(attempt)))
    MetadataForAttempt(attempt, withAttemptField)
  }

  private def attemptMetadataToIndexMetadata(index: ExecutionIndex, attemptMetadata: Iterable[MetadataForAttempt]) = {
    def addIndexProperty(value: JsObject) = JsObject(value.fields + (ShardKey -> JsNumber(index.getOrElse(-1))))
    val metadata = attemptMetadata.toList.sortBy(_.attempt) map { mdForAttempt => addIndexProperty(mdForAttempt.metadata) }
    MetadataForIndex(index.getOrElse(-1), metadata)
  }

  private def reduceWorkflowEvents(workflowEvents: Seq[MetadataEvent]): Seq[MetadataEvent] = {
    // This handles state specially so a sensible final value is returned irrespective of the order in which raw state
    // events were recorded in the journal.
    val (workflowStatusEvents, workflowNonStatusEvents) = workflowEvents partition(_.key.key == WorkflowMetadataKeys.Status)

    val ordering = implicitly[Ordering[WorkflowState]]
    // This orders by value in WorkflowState CRDT resolution, not necessarily the chronologically most recent state.
    val sortedStateEvents = workflowStatusEvents.filter(_.value.isDefined) sortWith { case (a, b) => ordering.gt(a.value.get.toWorkflowState, b.value.get.toWorkflowState) }
    workflowNonStatusEvents ++ sortedStateEvents.headOption.toList
  }

  private def parseWorkflowEventsToTimestampedJsValue(events: Seq[MetadataEvent]): JsObject = {
    // Partition if sequence of events in a pair of (Workflow level events, Call level events)
    val (workflowLevel, callLevel) = events partition { _.key.jobKey.isEmpty }
    val foldedWorkflowValues = eventsToIndexedJson(reduceWorkflowEvents(workflowLevel)).toJson.asJsObject

    val callsGroupedByFQN = callLevel groupBy { _.key.jobKey.get.callFqn }
    val callsGroupedByFQNAndIndex = callsGroupedByFQN mapValues { _ groupBy { _.key.jobKey.get.index } }
    val callsGroupedByFQNAndIndexAndAttempt = callsGroupedByFQNAndIndex mapValues { _ mapValues { _ groupBy { _.key.jobKey.get.attempt } } }

    val callsMap = callsGroupedByFQNAndIndexAndAttempt mapValues { eventsForIndex =>
      eventsForIndex mapValues { eventsForAttempt =>
        eventsForAttempt map Function.tupled(eventsToAttemptMetadata)
      } map { Function.tupled(attemptMetadataToIndexMetadata) }
    } mapValues { md => JsArray(md.toVector.sortBy(_.index) flatMap { _.metadata }) }

    val wrappedCalls = JsObject(Map(WorkflowMetadataKeys.Calls -> JsObject(callsMap)))

    JsObject(foldedWorkflowValues.fields ++ wrappedCalls.fields)
  }

  private def parseWorkflowEvents(events: Seq[MetadataEvent]): JsObject = parseWorkflowEventsToTimestampedJsValue(events)

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
    JsObject(Map(
      WorkflowMetadataKeys.Status -> JsString(status.toString),
      WorkflowMetadataKeys.Id -> JsString(workflowId.toString)
    ))
  }

  private def workflowMetadataResponse(workflowId: WorkflowId, eventsList: Seq[MetadataEvent]) = JsObject(MetadataBuilderActor.parseWorkflowEvents(eventsList).fields + ("id" -> JsString(workflowId.toString)))
}
