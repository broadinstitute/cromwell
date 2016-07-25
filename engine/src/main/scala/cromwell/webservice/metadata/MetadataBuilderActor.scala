package cromwell.webservice.metadata

import java.time.OffsetDateTime

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.core.Dispatcher.ApiDispatcher
import cromwell.core.ExecutionIndex.ExecutionIndex
import cromwell.core.{WorkflowId, WorkflowMetadataKeys, WorkflowState}
import cromwell.services.ServiceRegistryActor.ServiceRegistryFailure
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import cromwell.webservice.PerRequest.{RequestComplete, RequestCompleteWithHeaders}
import cromwell.webservice.metadata.IndexedJsonValue._
import cromwell.webservice.metadata.MetadataBuilderActor.{Idle, MetadataBuilderActorState, WaitingForMetadataService}
import cromwell.webservice.{APIResponse, WorkflowJsonSupport}
import org.slf4j.LoggerFactory
import spray.http.StatusCodes
import spray.httpx.SprayJsonSupport._
import spray.json._

import scala.collection.immutable.TreeMap
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import scalaz.std.list._
import scalaz.syntax.foldable._


object MetadataBuilderActor {
  sealed trait MetadataBuilderActorState
  case object Idle extends MetadataBuilderActorState
  case object WaitingForMetadataService extends MetadataBuilderActorState

  def props(serviceRegistryActor: ActorRef) = {
    Props(new MetadataBuilderActor(serviceRegistryActor)).withDispatcher(ApiDispatcher)
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
    // The `List` has a `Foldable` instance defined in scope, and because the `List`'s elements have a `Monoid` instance
    // defined in scope, `suml` can derive a sane `TimestampedJsValue` value even if the `List` of events is empty.
    events.toList map { e => keyValueToIndexedJson(e.key.key, e.value, e.offsetDateTime) } suml
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

  private def parseWorkflowEventsToTimestampedJsValue(events: Seq[MetadataEvent], includeCallsIfEmpty: Boolean): JsObject = {
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
    val callData = if (callsMap.isEmpty && !includeCallsIfEmpty) Nil else wrappedCalls.fields
    JsObject(foldedWorkflowValues.fields ++ callData)
  }

  private def parseWorkflowEvents(includeCallsIfEmpty: Boolean)(events: Seq[MetadataEvent]): JsObject = parseWorkflowEventsToTimestampedJsValue(events, includeCallsIfEmpty)

  /**
    * Parse a Seq of MetadataEvent into a full Json metadata response.
    */
  private def parse(events: Seq[MetadataEvent]): JsObject = {
    JsObject(events.groupBy(_.key.workflowId.toString) mapValues parseWorkflowEvents(includeCallsIfEmpty = true))
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
    case Event(failure: ServiceRegistryFailure, _) =>
      val response = APIResponse.fail(new RuntimeException("Can't find metadata service"))
      context.parent ! RequestComplete(StatusCodes.InternalServerError, response)
      allDone
    case Event(WorkflowQuerySuccess(uri, response, metadata), _) =>
      context.parent ! RequestCompleteWithHeaders(response, generateLinkHeaders(uri, metadata):_*)
      allDone
    case Event(failure: WorkflowQueryFailure, _) =>
      context.parent ! RequestComplete(StatusCodes.BadRequest, APIResponse.fail(failure.reason))
      allDone
    case Event(WorkflowOutputsResponse(id, events), _) =>
      // Add in an empty output event if there aren't already any output events.
      val hasOutputs = events collectFirst { case MetadataEvent(MetadataKey(_, _, key), _, _) if key == WorkflowMetadataKeys.Outputs => () } isDefined
      val updatedEvents = if (hasOutputs) events else MetadataEvent.empty(MetadataKey(id, None, WorkflowMetadataKeys.Outputs)) +: events
      context.parent ! RequestComplete(StatusCodes.OK, workflowMetadataResponse(id, updatedEvents, includeCallsIfEmpty = false))
      allDone
    case Event(LogsResponse(w, l), _) =>
      context.parent ! RequestComplete(StatusCodes.OK, workflowMetadataResponse(w, l, includeCallsIfEmpty = false))
      allDone
    case Event(failure: MetadataServiceFailure, _) =>
      context.parent ! RequestComplete(StatusCodes.InternalServerError, APIResponse.error(failure.reason))
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
        case MetadataQuery(w, _, _, _, _) => workflowMetadataResponse(w, eventsList)
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

  private def workflowMetadataResponse(workflowId: WorkflowId, eventsList: Seq[MetadataEvent], includeCallsIfEmpty: Boolean = true) = {
    JsObject(MetadataBuilderActor.parseWorkflowEvents(includeCallsIfEmpty)(eventsList).fields + ("id" -> JsString(workflowId.toString)))
  }
}
