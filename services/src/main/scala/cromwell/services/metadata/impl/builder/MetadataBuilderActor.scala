package cromwell.services.metadata.impl.builder

import java.time.OffsetDateTime
import java.util.concurrent.atomic.AtomicLong
import akka.actor.{ActorRef, LoggingFSM, PoisonPill, Props}
import cats.data.Validated.{Invalid, Valid}
import cats.implicits.catsSyntaxValidatedId
import common.collections.EnhancedCollections._
import common.validation.ErrorOr
import common.validation.ErrorOr._
import cromwell.services.metadata.impl.builder.MetadataComponent._
import cromwell.core.ExecutionIndex.ExecutionIndex
import cromwell.core._
import cromwell.services._
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import cromwell.services.metadata.impl.builder.MetadataBuilderActor._
import mouse.all._
import org.slf4j.LoggerFactory
import spray.json._

import java.time.temporal.ChronoUnit
import java.util.Currency
import scala.language.postfixOps

object MetadataBuilderActor {
  sealed trait MetadataBuilderActorState
  case object Idle extends MetadataBuilderActorState
  case object WaitingForMetadataService extends MetadataBuilderActorState
  case object WaitingForSubWorkflows extends MetadataBuilderActorState
  case object WaitingForSubWorkflowCost extends MetadataBuilderActorState

  sealed trait MetadataBuilderActorData

  case object IdleData extends MetadataBuilderActorData
  final case class HasWorkData(target: ActorRef, originalRequest: BuildMetadataJsonAction)
      extends MetadataBuilderActorData

  // Classes extending this trait are used to track state when the actor has launched child
  // actors to collect metadata for subworkflows. This class aggregates data as it comes in,
  // and builds the complete output when all subworkflow data is present. There's one child
  // class for plain metadata queries and one for cost queries.
  sealed trait EventsCollectorData extends MetadataBuilderActorData {
    val target: ActorRef
    val originalRequest: BuildMetadataJsonAction
    val originalQuery: MetadataQuery
    val originalEvents: Seq[MetadataEvent]
    val subWorkflowsMetadata: Map[String, JsValue]
    val waitFor: Int

    def isComplete = subWorkflowsMetadata.size == waitFor
  }

  final case class HasReceivedEventsData(target: ActorRef,
                                         originalRequest: BuildMetadataJsonAction,
                                         originalQuery: MetadataQuery,
                                         originalEvents: Seq[MetadataEvent],
                                         subWorkflowsMetadata: Map[String, JsValue],
                                         waitFor: Int
  ) extends EventsCollectorData {
    def withSubWorkflow(id: String, metadata: JsValue) =
      this.copy(subWorkflowsMetadata = subWorkflowsMetadata + ((id, metadata)))
  }

  final case class HasReceivedCostEventsData(target: ActorRef,
                                             originalRequest: BuildMetadataJsonAction,
                                             originalQuery: MetadataQuery,
                                             originalEvents: Seq[MetadataEvent],
                                             originalStatus: WorkflowState,
                                             subWorkflowsMetadata: Map[String, JsValue],
                                             waitFor: Int
  ) extends EventsCollectorData {
    def withSubWorkflow(id: String, metadata: JsValue) =
      this.copy(subWorkflowsMetadata = subWorkflowsMetadata + ((id, metadata)))
  }

  def props(readMetadataWorkerMaker: () => Props,
            metadataReadRowNumberSafetyThreshold: Int,
            isForSubworkflows: Boolean = false
  ) =
    Props(new MetadataBuilderActor(readMetadataWorkerMaker, metadataReadRowNumberSafetyThreshold, isForSubworkflows))

  val log = LoggerFactory.getLogger("MetadataBuilder")

  private val AttemptKey = "attempt"
  private val ShardKey = "shardIndex"

  final private val DefaultCurrency = Currency.getInstance("USD")

  /**
    * Metadata for a call attempt
    */
  private case class MetadataForAttempt(attempt: Int, metadata: JsObject)

  /**
    * Metadata objects of all attempts for one shard
    */
  private case class MetadataForIndex(index: Int, metadata: List[JsObject])

  /**
   * Extract the list of subworkflow ids from a list of metadata events
   */
  private def extractSubworkflowIds(events: Seq[MetadataEvent]): Seq[String] =
    events
      .collect {
        case MetadataEvent(key, value, _) if key.key.endsWith(CallMetadataKeys.SubWorkflowId) => value map { _.value }
      }
      .flatten
      .distinct

  private def extractFromJsAs[A: Manifest](js: JsObject, fieldName: String): Option[A] =
    js.fields.get(fieldName) match {
      case Some(a: A) => Some(a)
      case _ => None
    }

  private def eventsToAttemptMetadata(
    subWorkflowMetadata: Map[String, JsValue]
  )(attempt: Int, events: Seq[MetadataEvent]) = {
    val withAttemptField = JsObject(
      MetadataComponent(events, subWorkflowMetadata).toJson.asJsObject.fields + (AttemptKey -> JsNumber(attempt))
    )
    MetadataForAttempt(attempt, withAttemptField)
  }

  private def attemptMetadataToIndexMetadata(index: ExecutionIndex, attemptMetadata: Iterable[MetadataForAttempt]) = {
    def addIndexProperty(value: JsObject) = JsObject(value.fields + (ShardKey -> JsNumber(index.getOrElse(-1))))
    val metadata = attemptMetadata.toList.sortBy(_.attempt) map { mdForAttempt =>
      addIndexProperty(mdForAttempt.metadata)
    }
    MetadataForIndex(index.getOrElse(-1), metadata)
  }

  private def buildMetadataJson(events: Seq[MetadataEvent],
                                includeCallsIfEmpty: Boolean,
                                expandedValues: Map[String, JsValue]
  ): JsObject = {
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
     * This is guaranteed by the groupBy invariant and the fact that filter preserves the ordering. (See scala doc for groupBy and filter)
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
    val callsGroupedByFQNAndIndex = callsGroupedByFQN safeMapValues { _ groupBy { _.key.jobKey.get.index } }
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
    val callsGroupedByFQNAndIndexAndAttempt = callsGroupedByFQNAndIndex safeMapValues {
      _ safeMapValues { _ groupBy { _.key.jobKey.get.attempt } }
    }

    val eventsToAttemptFunction = Function.tupled(eventsToAttemptMetadata(expandedValues) _)
    val attemptToIndexFunction = (attemptMetadataToIndexMetadata _).tupled

    val callsMap = callsGroupedByFQNAndIndexAndAttempt safeMapValues {
      _ safeMapValues { _ map eventsToAttemptFunction } map attemptToIndexFunction
    } safeMapValues { md =>
      JsArray(md.toVector.sortBy(_.index) flatMap { _.metadata })
    }

    val wrappedCalls = JsObject(Map(WorkflowMetadataKeys.Calls -> JsObject(callsMap)))
    val callData = if (callsMap.isEmpty && !includeCallsIfEmpty) Nil else wrappedCalls.fields
    JsObject(workflowLevelJson.fields ++ callData)
  }

  private def parseWorkflowEvents(includeCallsIfEmpty: Boolean, expandedValues: Map[String, JsValue])(
    events: Seq[MetadataEvent]
  ): JsObject =
    buildMetadataJson(events, includeCallsIfEmpty, expandedValues)

  /**
    * Parse a Seq of MetadataEvent into a full Json metadata response.
    */
  private def parse(events: Seq[MetadataEvent], expandedValues: Map[String, JsValue]): JsObject =
    JsObject(
      events.groupBy(_.key.workflowId.toString) safeMapValues parseWorkflowEvents(includeCallsIfEmpty = true,
                                                                                  expandedValues
      )
    )

  val actorIdIterator = new AtomicLong(0)

  def uniqueActorName(workflowId: String): String =
    s"${getClass.getSimpleName}.${actorIdIterator.getAndIncrement()}-for-$workflowId"

  case class JobKeyAndGrouping(jobKey: MetadataJobKey, grouping: String)

  def makeSyntheticGroupedExecutionEvents(jobKeyAndGrouping: JobKeyAndGrouping,
                                          events: List[MetadataEvent]
  ): List[MetadataEvent] = {
    // The input list of `events` might be incoherent since some events that were logically generated may not (yet) have been
    // recorded in the database. This code is written defensively to check for a start date in the event list. If there isn't one,
    // just return the original list of events because we can't sanely construct a synthetic event.
    // The end date will correspond to the execution event with the largest start date.
    val startTimeEvents = events.filter(_.key.key.endsWith(":startTime"))
    if (startTimeEvents.isEmpty) {
      events
    } else {
      val oldestStartTimeKey = startTimeEvents.minBy(_.value collect { case MetadataValue(s, _) =>
        OffsetDateTime.parse(s).toEpochSecond
      } get)
      val newestStartTimeKey = startTimeEvents.maxBy(_.value collect { case MetadataValue(s, _) =>
        OffsetDateTime.parse(s).toEpochSecond
      } get)
      // Search for an end date event corresponding to the newest start date event. Use the same prefix as on the
      // newestStartDateEvent to search for it.
      val executionEventPrefix = newestStartTimeKey.key.key.takeWhile(_ != ':')
      val endTimeKey = s"$executionEventPrefix:endTime"
      val endTimeEvent = events.find(_.key.key == endTimeKey)

      val syntheticDescriptionKey = oldestStartTimeKey.key.copy(key = s"$executionEventPrefix:description")
      val syntheticStartTimeKey = oldestStartTimeKey.key.copy(key = s"$executionEventPrefix:startTime")

      // The start event may have had a different subscript than the end event so make sure the synthetic event that was
      // generated from these events uses a consistent subscript.
      List(
        oldestStartTimeKey.copy(key = syntheticStartTimeKey),
        oldestStartTimeKey.copy(key = syntheticDescriptionKey,
                                value = Option(MetadataValue(jobKeyAndGrouping.grouping, MetadataString))
        )
      ) ++ endTimeEvent.toList
    }
  }

  def groupEvents(events: Seq[MetadataEvent]): Seq[MetadataEvent] = {
    // Find all the MetadataKeys with executionEvents[x]:y names. We want to record the `x` and note whether we see
    // a `y` that has a `grouping` value.

    val (executionEvents, nonExecutionEvents) = events.partition(_.key.key.startsWith("executionEvents["))
    // Match any characters except a closing square brace.
    val executionEventKeyPatternRe = "[^]]+".r
    // Group execution events by their execution event keys.
    // The `get` inside the `groupBy` is safe because there must be a closing brace on the opening brace.
    val executionEventsByKeys = executionEvents groupBy { e =>
      (e.key.key.substring("executionEvents[".length) |> executionEventKeyPatternRe.findFirstIn).get
    }
    // "grouped" and "ungrouped" refer to the ":grouping" attribute that may be present in execution event metadata.
    val (groupedExecutionEvents, ungroupedExecutionEvents) = executionEventsByKeys partition { case (k, es) =>
      es.exists(_.key.key == s"executionEvents[$k]:grouping")
    }

    // (jobKey + grouping) => eeKey => executionEvents
    val groupedExecutionEventsByGrouping = groupedExecutionEvents groupBy {
      // The `get` on the `jobKey` is safe because execution events are all job-based. The outer `get` is safe because we already verified
      // the presence of this "grouping" attribute in this data set above.
      case (k, es) =>
        es.collectFirst {
          case e if e.key.key == s"executionEvents[$k]:grouping" =>
            JobKeyAndGrouping(e.key.jobKey.get, e.value.get.value)
        } get
    } map { case (jkg, m) =>
      jkg -> m.values.toList.flatten
    }

    // Tuplize the grouping function so it can operate on the List elements directly.
    val tupledGrouper = (makeSyntheticGroupedExecutionEvents _).tupled
    nonExecutionEvents ++ ungroupedExecutionEvents.values.toList.flatten ++ (groupedExecutionEventsByGrouping.toList flatMap tupledGrouper)
  }

  def processMetadataEvents(query: MetadataQuery,
                            eventsList: Seq[MetadataEvent],
                            expandedValues: Map[String, JsValue]
  ): JsObject =
    // Should we send back some message ? Or even fail the request instead ?
    if (eventsList.isEmpty) JsObject(Map.empty[String, JsValue])
    else {
      query match {
        case mq: MetadataQuery =>
          workflowMetadataResponse(mq.workflowId, eventsList, includeCallsIfEmpty = true, expandedValues)
        case _ => MetadataBuilderActor.parse(eventsList, expandedValues)
      }
    }

  def processStatusResponse(workflowId: WorkflowId, status: WorkflowState): JsObject =
    JsObject(
      Map(
        WorkflowMetadataKeys.Status -> JsString(status.toString),
        WorkflowMetadataKeys.Id -> JsString(workflowId.toString)
      )
    )

  def processLabelsResponse(workflowId: WorkflowId, labels: Map[String, String]): JsObject = {
    val jsLabels = labels map { case (k, v) => k -> JsString(v) }
    JsObject(
      Map(
        WorkflowMetadataKeys.Id -> JsString(workflowId.toString),
        WorkflowMetadataKeys.Labels -> JsObject(jsLabels)
      )
    )
  }

  def processOutputsResponse(id: WorkflowId, events: Seq[MetadataEvent]): JsObject = {
    // Add in an empty output event if there aren't already any output events.
    val hasOutputs = events exists { _.key.key.startsWith(WorkflowMetadataKeys.Outputs + ":") }
    val updatedEvents =
      if (hasOutputs) events else MetadataEvent.empty(MetadataKey(id, None, WorkflowMetadataKeys.Outputs)) +: events

    workflowMetadataResponse(id, updatedEvents, includeCallsIfEmpty = false, Map.empty)
  }

  def workflowMetadataResponse(workflowId: WorkflowId,
                               eventsList: Seq[MetadataEvent],
                               includeCallsIfEmpty: Boolean,
                               expandedValues: Map[String, JsValue]
  ): JsObject =
    JsObject(
      MetadataBuilderActor
        .parseWorkflowEvents(includeCallsIfEmpty, expandedValues)(eventsList)
        .fields + ("id" -> JsString(workflowId.toString))
    )

  /*
   * Attempt to the cost of a single call attempt from its metadata json, assuming the correct metadata
   * exists to allow that. We depend on vmStartTime and vmCostPerHour being present. We also use vmEndTime,
   * or fall back to the current time if it is absent.
   *
   * If the metadata needed to compute cost is missing, return 0. If the VM cost or any of the dates
   * can't be parsed, return an error. We will also return an error if we find a negative value for
   * vmCostPerHour - this indicates an error when generating the cost.
   */
  def computeCost(callName: String, jsVal: JsValue): ErrorOr[BigDecimal] =
    jsVal match {
      case jsCall: JsObject =>
        extractFromJsAs[JsString](jsCall, CallMetadataKeys.VmStartTime) map { startTimeVal =>
          val startTimeErrorOr = ErrorOr(OffsetDateTime.parse(startTimeVal.value))

          val endTimeErrorOr =
            extractFromJsAs[JsString](jsCall, CallMetadataKeys.VmEndTime)
              .map(v => ErrorOr(OffsetDateTime.parse(v.value)))
              .getOrElse(OffsetDateTime.now().validNel)

          val rawCostPerHour = extractFromJsAs[JsNumber](jsCall, CallMetadataKeys.VmCostPerHour).map(_.value)
          val costPerHourErrorOr = rawCostPerHour match {
            case Some(c: BigDecimal) if c >= 0 => c.validNel
            // A costPerHour < 0 indicates an error
            case Some(_: BigDecimal) =>
              val index =
                extractFromJsAs[JsNumber](jsCall, "shardIndex").map(_.value.toString).getOrElse("-1")
              val attempt =
                extractFromJsAs[JsNumber](jsCall, "attempt").map(_.value.toString).getOrElse("1")
              s"Couldn't find valid vmCostPerHour for ${List(callName, index, attempt).mkString(".")}".invalidNel
            case None => BigDecimal(0).validNel
          }

          for {
            start <- startTimeErrorOr
            end <- endTimeErrorOr
            costPerHour <- costPerHourErrorOr
            vmRuntimeInMillis = start.until(end, ChronoUnit.MILLIS).toDouble
          } yield (vmRuntimeInMillis / (1000 * 60 * 60)) * costPerHour
        } getOrElse (BigDecimal(0).validNel)
      case _ => BigDecimal(0).validNel
    }
}

class MetadataBuilderActor(readMetadataWorkerMaker: () => Props,
                           metadataReadRowNumberSafetyThreshold: Int,
                           isForSubworkflows: Boolean
) extends LoggingFSM[MetadataBuilderActorState, MetadataBuilderActorData]
    with DefaultJsonProtocol {

  import MetadataBuilderActor._

  startWith(Idle, IdleData)
  val tag = self.path.name

  when(Idle) { case Event(action: BuildMetadataJsonAction, IdleData) =>
    val readActor = context.actorOf(readMetadataWorkerMaker.apply())

    readActor ! action
    goto(WaitingForMetadataService) using HasWorkData(sender(), action)
  }

  private def allDone() = {
    context stop self
    stay()
  }

  when(WaitingForMetadataService) {
    case Event(StatusLookupResponse(w, status), HasWorkData(target, originalRequest)) =>
      target ! SuccessfulMetadataJsonResponse(originalRequest, processStatusResponse(w, status))
      allDone()
    case Event(LabelLookupResponse(w, labels), HasWorkData(target, originalRequest)) =>
      target ! SuccessfulMetadataJsonResponse(originalRequest, processLabelsResponse(w, labels))
      allDone()
    case Event(WorkflowOutputsResponse(id, events), HasWorkData(target, originalRequest)) =>
      target ! SuccessfulMetadataJsonResponse(originalRequest, processOutputsResponse(id, events))
      allDone()
    case Event(LogsResponse(w, l), HasWorkData(target, originalRequest)) =>
      target ! SuccessfulMetadataJsonResponse(originalRequest,
                                              workflowMetadataResponse(w, l, includeCallsIfEmpty = false, Map.empty)
      )
      allDone()
    case Event(CostResponse(w, s, m), HasWorkData(target, originalRequest)) =>
      processCostResponse(w, s, m, target, originalRequest)
    case Event(MetadataLookupResponse(query, metadata), HasWorkData(target, originalRequest)) =>
      processMetadataResponse(query, metadata, target, originalRequest)
    case Event(FetchFailedJobsMetadataLookupResponse(metadata), HasWorkData(target, originalRequest)) =>
      processFailedJobsMetadataResponse(metadata, target, originalRequest)
    case Event(MetadataLookupFailedTooLargeResponse(query, metadataSizeRows), HasWorkData(target, originalRequest)) =>
      val metadataTooLargeNumberOfRowsException =
        new MetadataTooLargeNumberOfRowsException(query.workflowId,
                                                  metadataSizeRows,
                                                  metadataReadRowNumberSafetyThreshold
        )
      target ! FailedMetadataJsonResponse(originalRequest, metadataTooLargeNumberOfRowsException)
      allDone()
    case Event(MetadataLookupFailedTimeoutResponse(query), HasWorkData(target, originalRequest)) =>
      val metadataTooLargeTimeoutException = new MetadataTooLargeTimeoutException(query.workflowId)
      target ! FailedMetadataJsonResponse(originalRequest, metadataTooLargeTimeoutException)
      allDone()
    case Event(FetchFailedTasksTimeoutResponse(workflowId), HasWorkData(target, originalRequest)) =>
      val metadataTooLargeTimeoutException = new MetadataTooLargeTimeoutException(workflowId)
      target ! FailedMetadataJsonResponse(originalRequest, metadataTooLargeTimeoutException)
      allDone()
    case Event(failure: MetadataServiceFailure, HasWorkData(target, originalRequest)) =>
      target ! FailedMetadataJsonResponse(originalRequest, failure.reason)
      allDone()
    case Event(response: RootAndSubworkflowLabelsLookupResponse, HasWorkData(target, _)) =>
      target ! response
      allDone()
  }

  when(WaitingForSubWorkflows) {
    case Event(mbr: MetadataJsonResponse, data: HasReceivedEventsData) =>
      processSubWorkflowMetadata(mbr, data)
    case Event(failure: MetadataServiceFailure, data: HasReceivedEventsData) =>
      data.target ! FailedMetadataJsonResponse(data.originalRequest, failure.reason)
      allDone()
  }

  when(WaitingForSubWorkflowCost) {
    case Event(mbr: MetadataJsonResponse, data: HasReceivedCostEventsData) =>
      processSubWorkflowCost(mbr, data)
    case Event(failure: MetadataServiceFailure, data: HasReceivedEventsData) =>
      data.target ! FailedMetadataJsonResponse(data.originalRequest, failure.reason)
      allDone()
  }

  whenUnhandled {
    case Event(message, IdleData) =>
      log.error(s"Received unexpected message $message in state $stateName with $IdleData")
      stay()
    case Event(message, HasWorkData(target, _)) =>
      log.error(s"Received unexpected message $message in state $stateName with target: $target")
      self ! PoisonPill
      stay()
    case Event(message, MetadataBuilderActor.HasReceivedEventsData(target, _, _, _, _, _)) =>
      log.error(s"Received unexpected message $message in state $stateName with target: $target")
      self ! PoisonPill
      stay()
  }

  def processSubWorkflowMetadata(metadataResponse: MetadataJsonResponse, data: HasReceivedEventsData) =
    metadataResponse match {
      case SuccessfulMetadataJsonResponse(GetMetadataAction(queryKey, _), js) =>
        val subId: WorkflowId = queryKey.workflowId
        val newData = data.withSubWorkflow(subId.toString, js)

        if (newData.isComplete) {
          buildAndStop(data.originalQuery,
                       data.originalEvents,
                       newData.subWorkflowsMetadata,
                       data.target,
                       data.originalRequest
          )
        } else {
          stay() using newData
        }
      case FailedMetadataJsonResponse(originalRequest, e) =>
        failAndDie(new RuntimeException(s"Failed to retrieve metadata for a sub workflow ($originalRequest)", e),
                   data.target,
                   data.originalRequest
        )

      case other =>
        val message =
          s"Programmer Error: MetadataBuilderActor expected subworkflow metadata response type but got ${other.getClass.getSimpleName}"
        log.error(message)
        failAndDie(new Exception(message), data.target, data.originalRequest)
    }

  def processSubWorkflowCost(metadataResponse: MetadataJsonResponse, data: HasReceivedCostEventsData) =
    metadataResponse match {
      case SuccessfulMetadataJsonResponse(GetCost(workflowId), js) =>
        val subId: WorkflowId = workflowId
        val newData = data.withSubWorkflow(subId.toString, js)

        if (newData.isComplete) {
          buildCostAndStop(
            data.originalQuery.workflowId,
            data.originalStatus,
            data.originalEvents,
            newData.subWorkflowsMetadata,
            data.target,
            data.originalRequest
          )
        } else {
          stay() using newData
        }
      case FailedMetadataJsonResponse(originalRequest, e) =>
        failAndDie(new RuntimeException(s"Failed to retrieve cost for a sub workflow ($originalRequest)", e),
                   data.target,
                   data.originalRequest
        )

      case other =>
        val message =
          s"Programmer Error: MetadataBuilderActor expected subworkflow metadata response type but got ${other.getClass.getSimpleName}"
        log.error(message)
        failAndDie(new Exception(message), data.target, data.originalRequest)
    }

  def failAndDie(reason: Throwable, target: ActorRef, originalRequest: BuildMetadataJsonAction) = {
    target ! FailedMetadataJsonResponse(originalRequest, reason)
    context stop self
    stay()
  }

  def buildAndStop(query: MetadataQuery,
                   eventsList: Seq[MetadataEvent],
                   expandedValues: Map[String, JsValue],
                   target: ActorRef,
                   originalRequest: BuildMetadataJsonAction
  ) = {
    val groupedEvents = groupEvents(eventsList)
    val res = processMetadataEvents(query, groupedEvents, expandedValues).fields
    target ! SuccessfulMetadataJsonResponse(originalRequest, JsObject(res))
    allDone()
  }

  def processMetadataResponse(query: MetadataQuery,
                              eventsList: Seq[MetadataEvent],
                              target: ActorRef,
                              originalRequest: BuildMetadataJsonAction
  ) =
    if (query.expandSubWorkflows) {
      // Scan events for sub workflow ids
      val subWorkflowIds = extractSubworkflowIds(eventsList)

      // If none is found just proceed to build metadata
      if (subWorkflowIds.isEmpty) buildAndStop(query, eventsList, Map.empty, target, originalRequest)
      else {
        // Otherwise spin up a metadata builder actor for each sub workflow
        subWorkflowIds foreach { subId =>
          val subMetadataBuilder = context.actorOf(MetadataBuilderActor.props(readMetadataWorkerMaker,
                                                                              metadataReadRowNumberSafetyThreshold,
                                                                              isForSubworkflows = true
                                                   ),
                                                   uniqueActorName(subId)
          )
          subMetadataBuilder ! GetMetadataAction(query.copy(workflowId = WorkflowId.fromString(subId)),
                                                 checkTotalMetadataRowNumberBeforeQuerying = false
          )
        }
        goto(WaitingForSubWorkflows) using HasReceivedEventsData(target,
                                                                 originalRequest,
                                                                 query,
                                                                 eventsList,
                                                                 Map.empty,
                                                                 subWorkflowIds.size
        )
      }
    } else {
      buildAndStop(query, eventsList, Map.empty, target, originalRequest)
    }

  def processCostResponse(id: WorkflowId,
                          status: WorkflowState,
                          metadataResponse: MetadataLookupResponse,
                          target: ActorRef,
                          originalRequest: BuildMetadataJsonAction
  ): State = {

    // Always expand subworkflows for cost
    val subWorkflowIds = extractSubworkflowIds(metadataResponse.eventList)

    if (subWorkflowIds.isEmpty)
      // If no subworkflows found, just build cost data
      buildCostAndStop(id, status, metadataResponse.eventList, Map.empty, target, originalRequest)
    else {
      // Otherwise spin up a metadata builder actor for each sub workflow
      subWorkflowIds foreach { subId =>
        val subMetadataBuilder = context.actorOf(MetadataBuilderActor.props(readMetadataWorkerMaker,
                                                                            metadataReadRowNumberSafetyThreshold,
                                                                            isForSubworkflows = true
                                                 ),
                                                 uniqueActorName(subId)
        )
        subMetadataBuilder ! GetCost(WorkflowId.fromString(subId))
      }
      goto(WaitingForSubWorkflowCost) using HasReceivedCostEventsData(target,
                                                                      originalRequest,
                                                                      metadataResponse.query,
                                                                      metadataResponse.eventList,
                                                                      status,
                                                                      Map.empty,
                                                                      subWorkflowIds.size
      )
    }
  }

  def buildCostAndStop(id: WorkflowId,
                       status: WorkflowState,
                       eventsList: Seq[MetadataEvent],
                       expandedValues: Map[String, JsValue],
                       target: ActorRef,
                       originalRequest: BuildMetadataJsonAction
  ): State = {

    val metadataEvents = MetadataBuilderActor.parse(groupEvents(eventsList), expandedValues)

    // Walk the structured metadata to attempt to compute cost for each call/shard/attempt in this workflow
    val callCostsWithErrors: List[ErrorOr[BigDecimal]] = for {
      wfObj <- extractFromJsAs[JsObject](metadataEvents, id.toString).toList
      callsObj <- extractFromJsAs[JsObject](wfObj, "calls").toList
      singleCallName <- callsObj.fields.keys.toList
      singleCallAttemptList <- extractFromJsAs[JsArray](callsObj, singleCallName).toList
      singleCallAttempt <- singleCallAttemptList.elements.toList
      singleCallAttemptCost = computeCost(singleCallName, singleCallAttempt)
    } yield singleCallAttemptCost

    val callCost: BigDecimal = callCostsWithErrors.collect { case (Valid(b)) => b }.sum
    val costErrors: Vector[JsString] =
      callCostsWithErrors
        .collect { case (Invalid(e)) => e.toList.mkString(", ") }
        .map(JsString(_))
        .toVector

    // Get the cost and error collection for any subworkflows that are part of this workflow
    val subworkflowCost = expandedValues.values.map {
      case o: JsObject =>
        extractFromJsAs[JsNumber](o, "cost").map(_.value).getOrElse(BigDecimal(0))
      case _ => BigDecimal(0)
    }.sum

    val subworkflowErrors = expandedValues.values.flatMap {
      case o: JsObject =>
        extractFromJsAs[JsArray](o, "errors").getOrElse(JsArray.empty).elements
      case _ => JsArray.empty.elements
    }

    val resp = JsObject(
      Map(
        WorkflowMetadataKeys.Id -> JsString(id.toString),
        WorkflowMetadataKeys.Status -> JsString(status.toString),
        "currency" -> JsString(DefaultCurrency.getCurrencyCode),
        "cost" -> JsNumber(callCost + subworkflowCost),
        "errors" -> JsArray(costErrors ++ subworkflowErrors)
      )
    )

    target ! SuccessfulMetadataJsonResponse(originalRequest, resp)
    allDone()
  }

  def processFailedJobsMetadataResponse(eventsList: Seq[MetadataEvent],
                                        target: ActorRef,
                                        originalRequest: BuildMetadataJsonAction
  ) = {
    val groupedEvents = groupEvents(eventsList)
    val res = MetadataBuilderActor.parse(groupedEvents, Map.empty)
    target ! SuccessfulMetadataJsonResponse(originalRequest, res)
    allDone()
  }
}
