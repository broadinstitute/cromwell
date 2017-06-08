package cromwell.webservice.metadata

import java.time.OffsetDateTime

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.webservice.metadata.MetadataComponent._
import cats.instances.list._
import cats.syntax.foldable._
import cromwell.core.Dispatcher.ApiDispatcher
import cromwell.core.ExecutionIndex.ExecutionIndex
import cromwell.core.{WorkflowId, WorkflowMetadataKeys, WorkflowState}
import cromwell.services.ServiceRegistryActor.ServiceRegistryFailure
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import cromwell.webservice.PerRequest.{RequestComplete, RequestCompleteWithHeaders}
import cromwell.webservice.metadata.MetadataBuilderActor.{Idle, MetadataBuilderActorData, MetadataBuilderActorState, WaitingForMetadataService, WaitingForSubWorkflows}
import cromwell.webservice.{APIResponse, PerRequestCreator, WorkflowJsonSupport}
import org.slf4j.LoggerFactory
import spray.http.{StatusCodes, Uri}
import spray.httpx.SprayJsonSupport._
import spray.json._

import scala.collection.immutable.TreeMap
import scala.language.postfixOps
import scala.util.{Failure, Random, Success, Try}


object MetadataBuilderActor {
  sealed trait MetadataBuilderActorState
  case object Idle extends MetadataBuilderActorState
  case object WaitingForMetadataService extends MetadataBuilderActorState
  case object WaitingForSubWorkflows extends MetadataBuilderActorState

  case class MetadataBuilderActorData(
                                      originalQuery: MetadataQuery,
                                      originalEvents: Seq[MetadataEvent],
                                      subWorkflowsMetadata: Map[String, JsValue],
                                      waitFor: Int
                                    ) {
    def withSubWorkflow(id: String, metadata: JsValue) = {
      this.copy(subWorkflowsMetadata = subWorkflowsMetadata + ((id, metadata)))
    }
    
    def isComplete = subWorkflowsMetadata.size == waitFor
  }
  
  def props(serviceRegistryActor: ActorRef) = {
    Props(new MetadataBuilderActor(serviceRegistryActor)).withDispatcher(ApiDispatcher)
  }

  val log = LoggerFactory.getLogger("MetadataBuilder")

  private val KeySeparator = MetadataKey.KeySeparator
  // Split on every unescaped KeySeparator
  private val KeySplitter = s"(?<!\\\\)$KeySeparator"
  private val bracketMatcher = """\[(\d*)\]""".r
  private val AttemptKey = "attempt"
  private val ShardKey = "shardIndex"

  def parseKeyChunk(chunk: String, innerValue: MetadataComponent): MetadataComponent = {
    chunk.indexOf('[') match {
      // If there's no bracket, it's an object. e.g.: "calls"
      case -1 => MetadataObject(Map(chunk -> innerValue))
      // If there's a bracket it's a named list. e.g.: "executionEvents[0][1]"
      case bracketIndex =>
        // Name: "executionEvents"
        val objectName = chunk.substring(0, bracketIndex)
        
        // Empty value means empty list
        if (innerValue == MetadataEmpty) MetadataObject(Map(objectName -> MetadataList(Map.empty)))
        else {
          // Brackets: "[0][1]"
          val brackets = chunk.substring(bracketIndex)
          // Indices as a list: List(0, 1)
          val listIndices = for {
            m <- bracketMatcher.findAllMatchIn(brackets)
            // It's possible for a bracket pair to be empty, in which case we just give it a random number
            asInt = if (m.group(1).isEmpty) Random.nextInt() else m.group(1).toInt
          } yield asInt
          // Fold into a MetadataList: MetadataList(0 -> MetadataList(1 -> innerValue))
          val metadataList = listIndices.toList.foldRight(innerValue)((index, acc) => MetadataList(TreeMap(index -> acc)))

          MetadataObject(Map(objectName -> metadataList))
        }
    }
  }

  private def metadataValueToComponent(value: Option[MetadataValue], customOrdering: Option[Ordering[MetadataPrimitive]]): MetadataComponent = {
    value map { someValue =>
      val coerced: Try[MetadataPrimitive] = someValue.valueType match {
        case MetadataInt => Try(MetadataPrimitive(JsNumber(someValue.value.toInt), customOrdering))
        case MetadataNumber => Try(MetadataPrimitive(JsNumber(someValue.value.toDouble), customOrdering))
        case MetadataBoolean => Try(MetadataPrimitive(JsBoolean(someValue.value.toBoolean), customOrdering))
        case MetadataString => Try(MetadataPrimitive(JsString(someValue.value), customOrdering))
      }

      coerced match {
        case Success(v) => v
        case Failure(e) =>
          log.warn(s"Failed to coerce ${someValue.value} to ${someValue.valueType}. Falling back to String.", e)
          MetadataPrimitive(JsString(someValue.value), customOrdering)
      }
    } getOrElse MetadataEmpty
  }

  def customOrdering(event: MetadataEvent): Option[Ordering[MetadataPrimitive]] = event match {
    case MetadataEvent(MetadataKey(_, Some(_), key), _, _) if key == CallMetadataKeys.ExecutionStatus => Option(MetadataPrimitive.ExecutionStatusOrdering)
    case MetadataEvent(MetadataKey(_, None, key), _, _) if key == WorkflowMetadataKeys.Status => Option(MetadataPrimitive.WorkflowStateOrdering)
    case _ => None
  }
  
  private def toMetadataComponent(subWorkflowMetadata: Map[String, JsValue])(event: MetadataEvent) = {
    import MetadataKey._
    
    lazy val originalKeyAndPrimitive = (event.key.key, metadataValueToComponent(event.value, customOrdering(event)))
    
    val keyAndPrimitive = if (event.key.key.endsWith(CallMetadataKeys.SubWorkflowId)) {
      (for {
        metadataValue <- event.value
        subWorkflowMetadata <- subWorkflowMetadata.get(metadataValue.value)
        keyWithSubWorkflowMetadata = event.key.key.replace(CallMetadataKeys.SubWorkflowId, CallMetadataKeys.SubWorkflowMetadata)
        subWorkflowComponent = MetadataPrimitive(subWorkflowMetadata)
      } yield (keyWithSubWorkflowMetadata, subWorkflowComponent)) getOrElse originalKeyAndPrimitive
    } else originalKeyAndPrimitive

    keyAndPrimitive._1.split(KeySplitter).map(_.unescapeMeta).foldRight(keyAndPrimitive._2)(parseKeyChunk)
  }

  /**
    * Metadata for a call attempt
    */
  private case class MetadataForAttempt(attempt: Int, metadata: JsObject)

  /**
    * Metadata objects of all attempts for one shard
    */
  private case class MetadataForIndex(index: Int, metadata: List[JsObject])

  implicit val dateTimeOrdering: Ordering[OffsetDateTime] = scala.Ordering.fromLessThan(_ isBefore _)

  /** Sort events by timestamp, transform them into TimestampedJsValues, and merge them together. */
  private def eventsToIndexedJson(events: Seq[MetadataEvent], subWorkflowMetadata: Map[String, JsValue]): MetadataComponent = {
    // The `List` has a `Foldable` instance defined in scope, and because the `List`'s elements have a `Monoid` instance
    // defined in scope, `combineAll` can derive a sane `TimestampedJsValue` value even if the `List` of events is empty.
    events.toList map toMetadataComponent(subWorkflowMetadata) combineAll
  }

  private def eventsToAttemptMetadata(subWorkflowMetadata: Map[String, JsValue])(attempt: Int, events: Seq[MetadataEvent]) = {
    val withAttemptField = JsObject(eventsToIndexedJson(events, subWorkflowMetadata).toJson.asJsObject.fields + (AttemptKey -> JsNumber(attempt)))
    MetadataForAttempt(attempt, withAttemptField)
  }

  private def attemptMetadataToIndexMetadata(index: ExecutionIndex, attemptMetadata: Iterable[MetadataForAttempt]) = {
    def addIndexProperty(value: JsObject) = JsObject(value.fields + (ShardKey -> JsNumber(index.getOrElse(-1))))
    val metadata = attemptMetadata.toList.sortBy(_.attempt) map { mdForAttempt => addIndexProperty(mdForAttempt.metadata) }
    MetadataForIndex(index.getOrElse(-1), metadata)
  }

  private def buildMetadataJson(events: Seq[MetadataEvent], includeCallsIfEmpty: Boolean, expandedValues: Map[String, JsValue]): JsObject = {
    // Partition events into workflow level and call level events
    val (workflowLevel, callLevel) = events partition { _.key.jobKey.isEmpty }
    val workflowLevelJson = eventsToIndexedJson(workflowLevel, Map.empty).toJson.asJsObject

    /*
     * Map(
     *    "fqn" -> Seq[Events],
     *    "fqn2" -> Seq[Events],
     *    ...
     * )
     * Note that groupBy will preserve the ordering of the events in the Seq, which means that as long as the DB sorts them by timestamp, we can always assume the last one is the newest one.
     * This is guaranteed by the groupBy invariant and the fact that filter preservers the ordering. (See scala doc for groupBy and filter)
     */
    val callsGroupedByFQN = callLevel groupBy { _.key.jobKey.get.callFqn }
    /*
     * Map(
     *    "fqn" -> Map( //Shard index
     *                    Option(0) -> Seq[Events],
     *                    Option(1) -> Seq[Events]
     *                    ...
     *             ),
     *    ...
     * )
     */
    val callsGroupedByFQNAndIndex = callsGroupedByFQN mapValues { _ groupBy { _.key.jobKey.get.index } }
    /*
     * Map(
     *    "fqn" -> Map(
     *                Option(0) -> Map( //Attempt
     *                                    1 -> Seq[Events],
     *                                    2 -> Seq[Events],
     *                                ...
     *                             ),   
     *                ...
     *             ),
     *    ...
     * )
     */
    val callsGroupedByFQNAndIndexAndAttempt = callsGroupedByFQNAndIndex mapValues { _ mapValues { _ groupBy { _.key.jobKey.get.attempt } } }
    
    val eventsToAttemptFunction = Function.tupled(eventsToAttemptMetadata(expandedValues) _)
    val attemptToIndexFunction = (attemptMetadataToIndexMetadata _).tupled

    val callsMap = callsGroupedByFQNAndIndexAndAttempt mapValues { _ mapValues { _ map eventsToAttemptFunction } map attemptToIndexFunction } mapValues { md => 
      JsArray(md.toVector.sortBy(_.index) flatMap { _.metadata })
    }

    val wrappedCalls = JsObject(Map(WorkflowMetadataKeys.Calls -> JsObject(callsMap)))
    val callData = if (callsMap.isEmpty && !includeCallsIfEmpty) Nil else wrappedCalls.fields
    JsObject(workflowLevelJson.fields ++ callData)
  }

  private def parseWorkflowEvents(includeCallsIfEmpty: Boolean, expandedValues: Map[String, JsValue])(events: Seq[MetadataEvent]): JsObject = {
    buildMetadataJson(events, includeCallsIfEmpty, expandedValues)
  }

  /**
    * Parse a Seq of MetadataEvent into a full Json metadata response.
    */
  private def parse(events: Seq[MetadataEvent], expandedValues: Map[String, JsValue]): JsObject = {
    JsObject(events.groupBy(_.key.workflowId.toString) mapValues parseWorkflowEvents(includeCallsIfEmpty = true, expandedValues))
  }
}

class MetadataBuilderActor(serviceRegistryActor: ActorRef) extends LoggingFSM[MetadataBuilderActorState, Option[MetadataBuilderActorData]]
  with DefaultJsonProtocol with WorkflowQueryPagination {

  import WorkflowJsonSupport._

  startWith(Idle, None)
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
    case Event(MetadataLookupResponse(query, metadata), None) =>
      processMetadataResponse(query, metadata)
    case Event(StatusLookupResponse(w, status), _) =>
      context.parent ! RequestComplete((StatusCodes.OK, processStatusResponse(w, status)))
      allDone
    case Event(failure: ServiceRegistryFailure, _) =>
      val response = APIResponse.fail(new RuntimeException("Can't find metadata service"))
      context.parent ! RequestComplete((StatusCodes.InternalServerError, response))
      allDone
    case Event(WorkflowQuerySuccess(uri: Uri, response, metadata), _) =>
      context.parent ! RequestCompleteWithHeaders(response, generateLinkHeaders(uri, metadata):_*)
      allDone
    case Event(failure: WorkflowQueryFailure, _) =>
      context.parent ! RequestComplete((StatusCodes.BadRequest, APIResponse.fail(failure.reason)))
      allDone
    case Event(WorkflowOutputsResponse(id, events), _) =>
      // Add in an empty output event if there aren't already any output events.
      val hasOutputs = events exists { _.key.key.startsWith(WorkflowMetadataKeys.Outputs + ":") }
      val updatedEvents = if (hasOutputs) events else MetadataEvent.empty(MetadataKey(id, None, WorkflowMetadataKeys.Outputs)) +: events
      context.parent ! RequestComplete((StatusCodes.OK, workflowMetadataResponse(id, updatedEvents, includeCallsIfEmpty = false, Map.empty)))
      allDone
    case Event(LogsResponse(w, l), _) =>
      context.parent ! RequestComplete((StatusCodes.OK, workflowMetadataResponse(w, l, includeCallsIfEmpty = false, Map.empty)))
      allDone
    case Event(failure: MetadataServiceFailure, _) =>
      context.parent ! RequestComplete((StatusCodes.InternalServerError, APIResponse.error(failure.reason)))
      allDone
    case Event(unexpectedMessage, stateData) =>
      val response = APIResponse.fail(new RuntimeException(s"MetadataBuilderActor $tag(WaitingForMetadataService, $stateData) got an unexpected message: $unexpectedMessage"))
      context.parent ! RequestComplete((StatusCodes.InternalServerError, response))
      context stop self
      stay()
  }
  
  when(WaitingForSubWorkflows) {
    case Event(RequestComplete(metadata), Some(data)) =>
      processSubWorkflowMetadata(metadata, data)
  }
  
  whenUnhandled {
    case Event(message, data) =>
      log.error(s"Received unexpected message $message in state $stateName with data $data")
      stay()
  }
  
  def processSubWorkflowMetadata(metadataResponse: Any, data: MetadataBuilderActorData) = {
    metadataResponse match {
      case (StatusCodes.OK, js: JsObject) =>
        js.fields.get(WorkflowMetadataKeys.Id) match {
          case Some(subId: JsString) =>
            val newData = data.withSubWorkflow(subId.value, js)
            
            if (newData.isComplete) {
              buildAndStop(data.originalQuery, data.originalEvents, newData.subWorkflowsMetadata)
            } else {
              stay() using Option(newData)
            }
          case _ => failAndDie(new RuntimeException("Received unexpected response while waiting for sub workflow metadata."))
        }
      case _ => failAndDie(new RuntimeException("Failed to retrieve metadata for a sub workflow."))
    }
  }
  
  def failAndDie(reason: Throwable) = {
    context.parent ! RequestComplete((StatusCodes.InternalServerError, APIResponse.error(reason)))
    context stop self
    stay()
  }
  
  def buildAndStop(query: MetadataQuery, eventsList: Seq[MetadataEvent], expandedValues: Map[String, JsValue]) = {
    context.parent ! RequestComplete((StatusCodes.OK, processMetadataEvents(query, eventsList, expandedValues)))
    allDone
  }
  
  def processMetadataResponse(query: MetadataQuery, eventsList: Seq[MetadataEvent]) = {
    if (query.expandSubWorkflows) {
      // Scan events for sub workflow ids
      val subWorkflowIds = eventsList.collect({
        case MetadataEvent(key, value, _) if key.key.endsWith(CallMetadataKeys.SubWorkflowId) => value map { _.value }
      }).flatten

      // If none is found just proceed to build metadata
      if (subWorkflowIds.isEmpty) buildAndStop(query, eventsList, Map.empty)
      else {
        // Otherwise spin up a metadata builder actor for each sub workflow
        subWorkflowIds foreach { subId =>
          val subMetadataBuilder = context.actorOf(MetadataBuilderActor.props(serviceRegistryActor), PerRequestCreator.endpointActorName)
          subMetadataBuilder ! GetMetadataQueryAction(query.copy(workflowId = WorkflowId.fromString(subId)))
        }
        goto(WaitingForSubWorkflows) using Option(MetadataBuilderActorData(query, eventsList, Map.empty, subWorkflowIds.size))
      }
    } else {
      buildAndStop(query, eventsList, Map.empty)
    }
  }

  def processMetadataEvents(query: MetadataQuery, eventsList: Seq[MetadataEvent], expandedValues: Map[String, JsValue]): JsObject = {
    // Should we send back some message ? Or even fail the request instead ?
    if (eventsList.isEmpty) JsObject(Map.empty[String, JsValue])
    else {
      query match {
        case MetadataQuery(w, _, _, _, _, _) => workflowMetadataResponse(w, eventsList, includeCallsIfEmpty = true, expandedValues)
        case _ => MetadataBuilderActor.parse(eventsList, expandedValues)
      }
    }
  }

  def processStatusResponse(workflowId: WorkflowId, status: WorkflowState): JsObject = {
    JsObject(Map(
      WorkflowMetadataKeys.Status -> JsString(status.toString),
      WorkflowMetadataKeys.Id -> JsString(workflowId.toString)
    ))
  }

  private def workflowMetadataResponse(workflowId: WorkflowId, eventsList: Seq[MetadataEvent], includeCallsIfEmpty: Boolean, expandedValues: Map[String, JsValue]) = {
    JsObject(MetadataBuilderActor.parseWorkflowEvents(includeCallsIfEmpty, expandedValues)(eventsList).fields + ("id" -> JsString(workflowId.toString)))
  }
}
