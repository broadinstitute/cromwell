package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{ActorRef, LoggingFSM, Props}
import cats.data.NonEmptyList
import cats.instances.list._
import cats.syntax.foldable._
import cromwell.core.WorkflowId
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffActor.{CallCacheDiffActorData, _}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffQueryParameter.CallCacheDiffQueryCall
import cromwell.services.metadata.CallMetadataKeys.CallCachingKeys
import cromwell.services.metadata.MetadataService.{GetMetadataQueryAction, MetadataLookupResponse, MetadataServiceKeyLookupFailed}
import cromwell.services.metadata._
import cromwell.webservice.APIResponse
import cromwell.webservice.PerRequest.RequestComplete
import cromwell.webservice.WorkflowJsonSupport._
import cromwell.webservice.metadata.MetadataComponent._
import cromwell.webservice.metadata._
import spray.http.StatusCodes
import spray.httpx.SprayJsonSupport._

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object CallCacheDiffActor {
  private val PlaceholderMissingHashValue = MetadataPrimitive(MetadataValue("Error: there is a hash entry for this key but the value is null !"))
  private val CallAAndBNotFoundException = new Exception("callA and callB have been run on a previous version of Cromwell on which this endpoint was not supported.")
  private val CallANotFoundException = new Exception("callA has been run on a previous version of Cromwell on which this endpoint was not supported.")
  private val CallBNotFoundException = new Exception("callB has been run on a previous version of Cromwell on which this endpoint was not supported.")

  sealed trait CallCacheDiffActorState
  case object Idle extends CallCacheDiffActorState
  case object WaitingForMetadata extends CallCacheDiffActorState

  sealed trait CallCacheDiffActorData
  case object CallCacheDiffNoData extends CallCacheDiffActorData
  case class CallCacheDiffWithRequest(queryA: MetadataQuery,
                                      queryB: MetadataQuery,
                                      responseA: Option[MetadataLookupResponse],
                                      responseB: Option[MetadataLookupResponse],
                                      replyTo: ActorRef
                                     ) extends CallCacheDiffActorData

  def props(serviceRegistryActor: ActorRef) = Props(new CallCacheDiffActor(serviceRegistryActor))
}

class CallCacheDiffActor(serviceRegistryActor: ActorRef) extends LoggingFSM[CallCacheDiffActorState, CallCacheDiffActorData] {
  startWith(Idle, CallCacheDiffNoData)

  when(Idle) {
    case Event(CallCacheDiffQueryParameter(callA, callB), CallCacheDiffNoData) =>
      val queryA = makeMetadataQuery(callA)
      val queryB = makeMetadataQuery(callB)
      serviceRegistryActor ! GetMetadataQueryAction(queryA)
      serviceRegistryActor ! GetMetadataQueryAction(queryB)
      goto(WaitingForMetadata) using CallCacheDiffWithRequest(queryA, queryB, None, None, sender())
  }

  when(WaitingForMetadata) {
    // First Response
    // Response A
    case Event(response: MetadataLookupResponse, data @ CallCacheDiffWithRequest(queryA, _, None, None, _)) if queryA == response.query =>
      stay() using data.copy(responseA = Option(response))
    // Response B
    case Event(response: MetadataLookupResponse, data @ CallCacheDiffWithRequest(_, queryB, None, None, _)) if queryB == response.query =>
      stay() using data.copy(responseB = Option(response))
    // Second Response
    // Response A
    case Event(response: MetadataLookupResponse, CallCacheDiffWithRequest(queryA, queryB, None, Some(responseB), replyTo)) if queryA == response.query =>
      buildDiffAndRespond(queryA, queryB, response, responseB, replyTo)
    // Response B
    case Event(response: MetadataLookupResponse, CallCacheDiffWithRequest(queryA, queryB, Some(responseA), None, replyTo)) if queryB == response.query =>
      buildDiffAndRespond(queryA, queryB, responseA, response, replyTo)
    case Event(MetadataServiceKeyLookupFailed(_, failure), data: CallCacheDiffWithRequest) =>
      data.replyTo ! RequestComplete((StatusCodes.InternalServerError, APIResponse.error(failure)))
      context stop self
      stay()
  }

  /**
    * Builds a response and sends it back as Json.
    * The response is structured in the following way
    * {
    *   "callA": {
    *     -- information about call A --
    *   },
    *   "callB": {
    *     -- information about call B --
    *   },
    *   "hashDifferential": [
    *     {
    *       "hash key": {
    *         "callA": -- hash value for call A, or null --,
    *         "callB": -- hash value for call B, or null --
    *       }
    *     },
    *     ...
    *   ]
    * }
    */
  private def buildDiffAndRespond(queryA: MetadataQuery,
                                  queryB: MetadataQuery,
                                  responseA: MetadataLookupResponse,
                                  responseB: MetadataLookupResponse,
                                  replyTo: ActorRef) = {

    val response = diffHashes(responseA.eventList, responseB.eventList) match {
      case Success(diff) => 
        val diffObject = MetadataObject(Map(
        "callA" -> makeCallInfo(queryA, responseA.eventList),
        "callB" -> makeCallInfo(queryB, responseB.eventList),
        "hashDifferential" -> diff
      ))
        
      RequestComplete((StatusCodes.OK, metadataComponentJsonWriter.write(diffObject).asJsObject))
      case Failure(f) => RequestComplete((StatusCodes.NotFound, APIResponse.error(f)))
    }
    
    replyTo ! response

    context stop self
    stay()
  }

  /**
    * Generates the "info" section of callA or callB
    */
  private def makeCallInfo(query: MetadataQuery, eventList: Seq[MetadataEvent]): MetadataComponent = {
    val callKey = MetadataObject(Map(
      "workflowId" -> MetadataPrimitive(MetadataValue(query.workflowId.toString)),
      "callFqn" -> MetadataPrimitive(MetadataValue(query.jobKey.get.callFqn)),
      "jobIndex" -> MetadataPrimitive(MetadataValue(query.jobKey.get.index.getOrElse(-1)))
    ))

    val allowResultReuse = attributeToComponent(eventList, { _ == CallCachingKeys.AllowReuseMetadataKey }, { _ => "allowResultReuse" })
    val executionStatus = attributeToComponent(eventList, { _ == CallMetadataKeys.ExecutionStatus })

    List(callKey, allowResultReuse, executionStatus) combineAll
  }

  /**
    * Collects events from the list for which the keys verify the keyFilter predicate
    * and apply keyModifier to the event's key
    */
  private def collectEvents(events: Seq[MetadataEvent],
                            keyFilter: (String  => Boolean),
                            keyModifier: (String => String)) = events collect {
    case event @ MetadataEvent(metadataKey @ MetadataKey(_, _, key), _, _) if keyFilter(key) =>
      event.copy(key = metadataKey.copy(key = keyModifier(key)))
  }

  /**
    * Given a list of events, a keyFilter and a keyModifier, returns the associated MetadataComponent.
    * Ensures that events are properly aggregated together (CRDTs and latest timestamp rule)
    */
  private def attributeToComponent(events: Seq[MetadataEvent], keyFilter: (String  => Boolean), keyModifier: (String => String) = identity[String]) = {
    MetadataComponent(collectEvents(events, keyFilter, keyModifier))
  }

  /**
    * Makes a diff object out of a key and a pair of values.
    * Values are Option[Option[MetadataValue]] for the following reason:
    *
    * The outer option represents whether or not this key had a corresponding hash metadata entry for the given call
    * If the above is true, the inner value is the metadata value for this entry, which is nullable, hence an Option.
    * The first outer option will determine whether the resulting json value will be null (no hash entry for this key),
    * or the actual value.
    * If the metadata value (inner option) happens to be None, it's an error, as we don't expect to publish null hash values.
    * In that case we replace it with the placeholderMissingHashValue.
    */
  private def makeHashDiffObject(key: String, valueA: Option[Option[MetadataValue]], valueB: Option[Option[MetadataValue]]) = {
    def makeFinalValue(value: Option[Option[MetadataValue]]) = value match {
      case Some(Some(metadataValue)) => MetadataPrimitive(metadataValue)
      case Some(None) => PlaceholderMissingHashValue
      case None => MetadataNull
    }

    MetadataObject(key.trim ->
      MetadataObject(
        "callA" -> makeFinalValue(valueA),
        "callB" -> makeFinalValue(valueB)
      )
    )
  }

  /**
    * Creates the hash differential between 2 list of events
    */
  private def diffHashes(eventsA: Seq[MetadataEvent], eventsB: Seq[MetadataEvent]): Try[MetadataComponent] = {
    val hashesKey = CallCachingKeys.HashesKey + MetadataKey.KeySeparator
    // Collect hashes events and map their key to only keep the meaningful part of the key
    // Then map the result to get a Map of hashKey -> Option[MetadataValue]. This will allow for fast lookup when
    // comparing the 2 hash sets.
    // Note that it's an Option[MetadataValue] because metadata values can be null, although for this particular
    // case we don't expect it to be (we should never publish a hash metadata event with a null value)
    // If that happens we will place a placeholder value in place of the hash to signify of the unexpected absence of it
    def collectHashes(events: Seq[MetadataEvent]) = {
      collectEvents(events, { _.startsWith(hashesKey) }, { _.stripPrefix(hashesKey) })  map {
        case MetadataEvent(MetadataKey(_, _, keyA), valueA, _) => keyA -> valueA
      } toMap
    }

    val hashesA: Map[String, Option[MetadataValue]] = collectHashes(eventsA)
    val hashesB: Map[String, Option[MetadataValue]] = collectHashes(eventsB)

    (hashesA.isEmpty, hashesB.isEmpty) match {
      case (true, true) => Failure(CallAAndBNotFoundException)
      case (true, false) => Failure(CallANotFoundException)
      case (false, true) => Failure(CallBNotFoundException)
      case (false, false) => Success(diffHashEvents(hashesA, hashesB))
    }

  }
  
  private def diffHashEvents(hashesA: Map[String, Option[MetadataValue]], hashesB: Map[String, Option[MetadataValue]]) = {
    val hashesUniqueToB: Map[String, Option[MetadataValue]] = hashesB.filterNot({ case (k, _) => hashesA.keySet.contains(k) })

    val hashDiff: List[MetadataComponent] = {
      // Start with all hashes in A
      hashesA
        // Try to find the corresponding pair in B.
        // We end up with a 
        // List[(Option[String, Option[MetadataValue], Option[String, Option[MetadataValue])]
        //                ^               ^                     ^               ^
        //             hashKey        hashValue              hashKey         hashValue
        //               for             for                   for              for
        //                A               A                     B                B
        //      |____________________________________| |___________________________________|
        //                  hashPair for A                        hashPair for B
        // 
        // HashPairs are Some or None depending on whether or not they have a metadata entry for the corresponding hashKey
        // At this stage we only have Some(hashPair) for A, and either Some(hashPair) or None for B depending on if we found it in hashesB
        .map({
        hashPairA => Option(hashPairA) -> hashesB.find(_._1 == hashPairA._1)
      })
        // Add the missing hashes that are in B but not in A. The left hashPair is therefore None
        .++(hashesUniqueToB.map(None -> Option(_)))
        .collect({
          // Both have a value but they're different. We can assume the keys are the same (if we did our job right until here)
          case (Some((keyA, valueA)), Some((_, valueB))) if valueA != valueB =>
            makeHashDiffObject(keyA, Option(valueA), Option(valueB))
          // Key is in A but not in B
          case (Some((keyA, valueA)), None) =>
            makeHashDiffObject(keyA, Option(valueA), None)
          // Key is in B but not in A
          case (None, Some((keyB, valueB))) =>
            makeHashDiffObject(keyB, None, Option(valueB))
        })
        .toList
    }

    MetadataList(hashDiff)
  }

  /**
    * Create a Metadata query from a CallCacheDiffQueryCall
    */
  private def makeMetadataQuery(call: CallCacheDiffQueryCall) = MetadataQuery(
    WorkflowId.fromString(call.workflowId),
    // jobAttempt None will return keys for all attempts
    Option(MetadataQueryJobKey(call.callFqn, call.jobIndex, None)),
    None,
    Option(NonEmptyList.of("callCaching", "executionStatus")),
    None,
    expandSubWorkflows = false
  )
}
