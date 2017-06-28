package cromwell.webservice.metadata

import java.util.UUID

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.webservice.metadata.MetadataComponent._
import cromwell.core.Dispatcher.ApiDispatcher
import cromwell.core.ExecutionIndex.ExecutionIndex
import cromwell.core.{WorkflowId, WorkflowMetadataKeys, WorkflowState}
import cromwell.services.ServiceRegistryActor.ServiceRegistryFailure
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import cromwell.webservice.metadata.MetadataBuilderActor._
import org.slf4j.LoggerFactory
import spray.json._


object MetadataBuilderActor {
  sealed abstract class MetadataBuilderActorResponse
  case class BuiltMetadataResponse(response: JsObject) extends MetadataBuilderActorResponse
  case class FailedMetadataResponse(reason: Throwable) extends MetadataBuilderActorResponse

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

  private val AttemptKey = "attempt"
  private val ShardKey = "shardIndex"

  /**
    * Metadata for a call attempt
    */
  private case class MetadataForAttempt(attempt: Int, metadata: JsObject)

  /**
    * Metadata objects of all attempts for one shard
    */
  private case class MetadataForIndex(index: Int, metadata: List[JsObject])

  private def eventsToAttemptMetadata(subWorkflowMetadata: Map[String, JsValue])(attempt: Int, events: Seq[MetadataEvent]) = {
    val withAttemptField = JsObject(MetadataComponent(events, subWorkflowMetadata).toJson.asJsObject.fields + (AttemptKey -> JsNumber(attempt)))
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
    val workflowLevelJson = MetadataComponent(workflowLevel, Map.empty).toJson.asJsObject

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

  def uniqueActorName: String = List("MetadataBuilderActor", UUID.randomUUID()).mkString("-")
}

class MetadataBuilderActor(serviceRegistryActor: ActorRef) extends LoggingFSM[MetadataBuilderActorState, Option[MetadataBuilderActorData]]
  with DefaultJsonProtocol {
  import MetadataBuilderActor._

  private var target: ActorRef = ActorRef.noSender

  startWith(Idle, None)
  val tag = self.path.name

  when(Idle) {
    case Event(action: MetadataServiceAction, _) =>
      target = sender()
      serviceRegistryActor ! action
      goto(WaitingForMetadataService)
  }

  private def allDone = {
    context stop self
    stay()
  }

  when(WaitingForMetadataService) {
    case Event(StatusLookupResponse(w, status), _) =>
      target ! BuiltMetadataResponse(processStatusResponse(w, status))
      allDone
    case Event(WorkflowOutputsResponse(id, events), _) =>
      // Add in an empty output event if there aren't already any output events.
      val hasOutputs = events exists { _.key.key.startsWith(WorkflowMetadataKeys.Outputs + ":") }
      val updatedEvents = if (hasOutputs) events else MetadataEvent.empty(MetadataKey(id, None, WorkflowMetadataKeys.Outputs)) +: events
      target ! BuiltMetadataResponse(workflowMetadataResponse(id, updatedEvents, includeCallsIfEmpty = false, Map.empty))
      allDone
    case Event(LogsResponse(w, l), _) =>
      target ! BuiltMetadataResponse(workflowMetadataResponse(w, l, includeCallsIfEmpty = false, Map.empty))
      allDone
    case Event(MetadataLookupResponse(query, metadata), None) => processMetadataResponse(query, metadata)
    case Event(failure: ServiceRegistryFailure, _) =>
      target ! FailedMetadataResponse(new RuntimeException("Can't find metadata service"))
      allDone
    case Event(failure: MetadataServiceFailure, _) =>
      target ! FailedMetadataResponse(failure.reason)
      allDone
    case Event(unexpectedMessage, stateData) =>
      target ! FailedMetadataResponse(new RuntimeException(s"MetadataBuilderActor $tag(WaitingForMetadataService, $stateData) got an unexpected message: $unexpectedMessage"))
      context stop self
      stay()
  }

  when(WaitingForSubWorkflows) {
    case Event(mbr: MetadataBuilderActorResponse, Some(data)) =>
      processSubWorkflowMetadata(mbr, data)
  }

  whenUnhandled {
    case Event(message, data) =>
      log.error(s"Received unexpected message $message in state $stateName with data $data")
      stay()
  }

  def processSubWorkflowMetadata(metadataResponse: MetadataBuilderActorResponse, data: MetadataBuilderActorData) = {
    metadataResponse match {
      case BuiltMetadataResponse(js) =>
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
      case FailedMetadataResponse(e) => failAndDie(new RuntimeException("Failed to retrieve metadata for a sub workflow.", e))
    }
  }

  def failAndDie(reason: Throwable) = {
    target ! FailedMetadataResponse(reason)
    context stop self
    stay()
  }

  def buildAndStop(query: MetadataQuery, eventsList: Seq[MetadataEvent], expandedValues: Map[String, JsValue]) = {
    target ! BuiltMetadataResponse(processMetadataEvents(query, eventsList, expandedValues))
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
          val subMetadataBuilder = context.actorOf(MetadataBuilderActor.props(serviceRegistryActor), uniqueActorName)
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

  private def workflowMetadataResponse(workflowId: WorkflowId,
                                       eventsList: Seq[MetadataEvent],
                                       includeCallsIfEmpty: Boolean,
                                       expandedValues: Map[String, JsValue]): JsObject = {
    JsObject(MetadataBuilderActor.parseWorkflowEvents(includeCallsIfEmpty, expandedValues)(eventsList).fields + ("id" -> JsString(workflowId.toString)))
  }
}
