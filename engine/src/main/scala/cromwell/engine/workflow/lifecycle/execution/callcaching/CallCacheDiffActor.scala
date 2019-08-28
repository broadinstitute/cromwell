package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{ActorRef, LoggingFSM, Props}
import cats.data.NonEmptyList
import cats.syntax.validated._
import cats.syntax.apply._
import cats.syntax.traverse._
import cats.instances.list._
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import common.validation.Validation._
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffActor.{CallCacheDiffActorData, _}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffQueryParameter.CallCacheDiffQueryCall
import cromwell.services.metadata.MetadataService.{GetMetadataAction, MetadataServiceKeyLookupFailed}
import cromwell.services.metadata._
import cromwell.services.metadata.impl.builder.MetadataBuilderActor.BuiltMetadataResponse
import spray.json.{JsArray, JsObject, JsString, JsValue}

class CallCacheDiffActor(serviceRegistryActor: ActorRef) extends LoggingFSM[CallCacheDiffActorState, CallCacheDiffActorData] {
  startWith(Idle, CallCacheDiffNoData)

  when(Idle) {
    case Event(CallCacheDiffQueryParameter(callA, callB), CallCacheDiffNoData) =>
      val queryA = makeMetadataQuery(callA)
      val queryB = makeMetadataQuery(callB)
      serviceRegistryActor ! GetMetadataAction(queryA)
      serviceRegistryActor ! GetMetadataAction(queryB)
      goto(WaitingForMetadata) using CallCacheDiffWithRequest(queryA, queryB, None, None, sender())
  }

  when(WaitingForMetadata) {
    // First Response
    // Response A
    case Event(BuiltMetadataResponse(GetMetadataAction(originalQuery), responseJson), data@CallCacheDiffWithRequest(queryA, _, None, None, _)) if queryA == originalQuery =>
      stay() using data.copy(responseA = Option(WorkflowMetadataJson(responseJson)))
    // Response B
    case Event(BuiltMetadataResponse(GetMetadataAction(originalQuery), responseJson), data@CallCacheDiffWithRequest(_, queryB, None, None, _)) if queryB == originalQuery =>
      stay() using data.copy(responseB = Option(WorkflowMetadataJson(responseJson)))
    // Second Response
    // Response A
    case Event(BuiltMetadataResponse(GetMetadataAction(originalQuery), responseJson), CallCacheDiffWithRequest(queryA, queryB, None, Some(responseB), replyTo)) if queryA == originalQuery =>
      buildDiffAndRespond(queryA, queryB, WorkflowMetadataJson(responseJson), responseB, replyTo)
    // Response B
    case Event(BuiltMetadataResponse(GetMetadataAction(originalQuery), responseJson), CallCacheDiffWithRequest(queryA, queryB, Some(responseA), None, replyTo)) if queryB == originalQuery =>
      buildDiffAndRespond(queryA, queryB, responseA, WorkflowMetadataJson(responseJson), replyTo)
    case Event(MetadataServiceKeyLookupFailed(_, failure), data: CallCacheDiffWithRequest) =>
      data.replyTo ! FailedCallCacheDiffResponse(failure)
      context stop self
      stay()
  }

  whenUnhandled {
    case Event(oops, oopsData) =>
      log.error(s"Programmer Error: Unexpected event received by ${this.getClass.getSimpleName}: $oops / $oopsData (in state $stateName)")
      stay()

  }

  /**
    * Builds a response and sends it back as Json.
    * The response is structured in the following way
    * {
    * "callA": {
    * -- information about call A --
    * },
    * "callB": {
    * -- information about call B --
    * },
    * "hashDifferential": [
    * {
    * "hash key": {
    * "callA": -- hash value for call A, or null --,
    * "callB": -- hash value for call B, or null --
    * }
    * },
    * ...
    * ]
    * }
    */
  private def buildDiffAndRespond(queryA: MetadataQuery,
                                  queryB: MetadataQuery,
                                  responseA: WorkflowMetadataJson,
                                  responseB: WorkflowMetadataJson,
                                  replyTo: ActorRef) = {

    //    lazy val buildResponse = {
    //      diffHashes(responseA.responseJson, responseB.responseJson) match {
    //        case Success(diff) =>
    //          val diffObject = MetadataObject(Map(
    //            "callA" -> makeCallInfo(queryA, responseA.responseJson),
    //            "callB" -> makeCallInfo(queryB, responseB.responseJson),
    //            "hashDifferential" -> diff
    //          ))
    //
    //          BuiltCallCacheDiffResponse(metadataComponentJsonWriter.write(diffObject).asJsObject)
    //        case Failure(f) => FailedCallCacheDiffResponse(f)
    //      }
    //    }
    //
    //    val response = checkCallsExistence(queryA, queryB, responseA, responseB) match {
    //      case Some(msg) => FailedCallCacheDiffResponse(CachedCallNotFoundException(msg))
    //      case None => buildResponse
    //    }

    val callACachingMetadata = extractCallCachingMetadata(queryA, responseA)
    val callBCachingMetadata = extractCallCachingMetadata(queryB, responseB)

    (callACachingMetadata, callBCachingMetadata) mapN { case (callA, callB) =>
      println(extractHashes(callA))
      println(extractHashes(callB))
    }

    val response = FailedCallCacheDiffResponse(CachedCallNotFoundException("Not implemented"))

    replyTo ! response

    context stop self
    stay()
  }

  //  /**
  //    * Returns an error message if one or both of the calls are not found, or None if it does
  //    */
  //  private def checkCallsExistence(queryA: MetadataQuery,
  //                                  queryB: MetadataQuery,
  //                                  responseA: BuiltMetadataResponse,
  //                                  responseB: BuiltMetadataResponse): Option[String] = {
  //    import cromwell.core.ExecutionIndex._
  //
  //    def makeTag(query: MetadataQuery) = {
  //      s"${query.workflowId}:${query.jobKey.get.callFqn}:${query.jobKey.get.index.fromIndex}"
  //    }
  //
  //    def makeNotFoundMessage(queries: NonEmptyList[MetadataQuery]) = {
  //      val plural = if (queries.tail.nonEmpty) "s" else ""
  //      s"Cannot find call$plural ${queries.map(makeTag).toList.mkString(", ")}"
  //    }
  //
  //    (responseA.eventList, responseB.eventList) match {
  //      case (a, b) if a.isEmpty && b.isEmpty => Option(makeNotFoundMessage(NonEmptyList.of(queryA, queryB)))
  //      case (a, _) if a.isEmpty => Option(makeNotFoundMessage(NonEmptyList.of(queryA)))
  //      case (_, b) if b.isEmpty => Option(makeNotFoundMessage(NonEmptyList.of(queryB)))
  //      case _ => None
  //    }
  //  }

  //  /**
  //    * Generates the "info" section of callA or callB
  //    */
  //  private def makeCallInfo(query: MetadataQuery, eventList: JsObject): MetadataComponent = {
  //    val callKey = MetadataObject(Map(
  //      "workflowId" -> MetadataPrimitive(MetadataValue(query.workflowId.toString)),
  //      "callFqn" -> MetadataPrimitive(MetadataValue(query.jobKey.get.callFqn)),
  //      "jobIndex" -> MetadataPrimitive(MetadataValue(query.jobKey.get.index.getOrElse(-1)))
  //    ))
  //
  //    val allowResultReuse = attributeToComponent(eventList, { _ == CallCachingKeys.AllowReuseMetadataKey }, { _ => "allowResultReuse" })
  //    val executionStatus = attributeToComponent(eventList, { _ == CallMetadataKeys.ExecutionStatus })
  //
  //    List(callKey, allowResultReuse, executionStatus) combineAll
  //  }
  //
  //  /**
  //    * Collects events from the list for which the keys verify the keyFilter predicate
  //    * and apply keyModifier to the event's key
  //    */
  //  private def collectEvents(events: Seq[MetadataEvent],
  //                            keyFilter: (String  => Boolean),
  //                            keyModifier: (String => String)) = events collect {
  //    case event @ MetadataEvent(metadataKey @ MetadataKey(_, _, key), _, _) if keyFilter(key) =>
  //      event.copy(key = metadataKey.copy(key = keyModifier(key)))
  //  }

  //  /**
  //    * Given a list of events, a keyFilter and a keyModifier, returns the associated MetadataComponent.
  //    * Ensures that events are properly aggregated together (CRDTs and latest timestamp rule)
  //    */
  //  private def attributeToComponent(events: Seq[MetadataEvent], keyFilter: (String  => Boolean), keyModifier: (String => String) = identity[String]) = {
  //    MetadataComponent(collectEvents(events, keyFilter, keyModifier))
  //  }

  //  /**
  //    * Makes a diff object out of a key and a pair of values.
  //    * Values are Option[Option[MetadataValue]] for the following reason:
  //    *
  //    * The outer option represents whether or not this key had a corresponding hash metadata entry for the given call
  //    * If the above is true, the inner value is the metadata value for this entry, which is nullable, hence an Option.
  //    * The first outer option will determine whether the resulting json value will be null (no hash entry for this key),
  //    * or the actual value.
  //    * If the metadata value (inner option) happens to be None, it's an error, as we don't expect to publish null hash values.
  //    * In that case we replace it with the placeholderMissingHashValue.
  //    */
  //  private def makeHashDiffObject(key: String, valueA: Option[Option[MetadataValue]], valueB: Option[Option[MetadataValue]]) = {
  //    def makeFinalValue(value: Option[Option[MetadataValue]]) = value match {
  //      case Some(Some(metadataValue)) => MetadataPrimitive(metadataValue)
  //      case Some(None) => PlaceholderMissingHashValue
  //      case None => MetadataNullComponent
  //    }
  //
  //    MetadataObject(
  //      "hashKey" -> MetadataPrimitive(MetadataValue(key.trim, MetadataString)),
  //      "callA" -> makeFinalValue(valueA),
  //      "callB" -> makeFinalValue(valueB)
  //    )
  //  }

}


object CallCacheDiffActor {
  //  private val PlaceholderMissingHashValue = MetadataPrimitive(MetadataValue("Error: there is a hash entry for this key but the value is null !"))

  final case class CachedCallNotFoundException(message: String) extends Exception {
    override def getMessage = message
  }

  // Exceptions when calls exist but have no hashes in their metadata, indicating they were run pre-28
  //  private val HashesForCallAAndBNotFoundException = new Exception("callA and callB have not finished yet, or were run on a previous version of Cromwell on which this endpoint was not supported.")
  //  private val HashesForCallANotFoundException = new Exception("callA has not finished yet, or was run on a previous version of Cromwell on which this endpoint was not supported.")
  //  private val HashesForCallBNotFoundException = new Exception("callB has not finished yet, or was run on a previous version of Cromwell on which this endpoint was not supported.")

  sealed trait CallCacheDiffActorState
  case object Idle extends CallCacheDiffActorState
  case object WaitingForMetadata extends CallCacheDiffActorState

  sealed trait CallCacheDiffActorData
  case object CallCacheDiffNoData extends CallCacheDiffActorData
  case class CallCacheDiffWithRequest(queryA: MetadataQuery,
                                      queryB: MetadataQuery,
                                      responseA: Option[WorkflowMetadataJson],
                                      responseB: Option[WorkflowMetadataJson],
                                      replyTo: ActorRef
                                     ) extends CallCacheDiffActorData

  sealed abstract class CallCacheDiffActorResponse
  case class BuiltCallCacheDiffResponse(response: JsObject) extends CallCacheDiffActorResponse
  case class FailedCallCacheDiffResponse(reason: Throwable) extends CallCacheDiffActorResponse


  def props(serviceRegistryActor: ActorRef) = Props(new CallCacheDiffActor(serviceRegistryActor)).withDispatcher(EngineDispatcher)

  /**
    * Create a Metadata query from a CallCacheDiffQueryCall
    */
  def makeMetadataQuery(call: CallCacheDiffQueryCall) = MetadataQuery(
    workflowId = call.workflowId,
    // jobAttempt None will return keys for all attempts
    jobKey = Option(MetadataQueryJobKey(call.callFqn, call.jobIndex, None)),
    key = None,
    includeKeysOption = Option(NonEmptyList.of("callCaching", "executionStatus")),
    excludeKeysOption = Option(NonEmptyList.of("callCaching:hitFailures")),
    expandSubWorkflows = false
  )

  // These simple case classes are just to help apply a little type safety to input and output types:
  final case class WorkflowMetadataJson(value: JsObject) extends AnyVal
  final case class CallCachingMetadataJson(value: JsObject) extends AnyVal

  /*
   * Takes in the JsObject returned from a metadata query and filters out only the appropriate call's callCaching section
   */
  def extractCallCachingMetadata(query: MetadataQuery, response: WorkflowMetadataJson): ErrorOr[CallCachingMetadataJson] = {

    for {
      // Sanity Checks:
      _ <- response.value.checkFieldValue("id", s""""${query.workflowId}"""")
      jobKey <- query.jobKey.toErrorOr("Call is required in call cache diff query")

      // Unpack the JSON:
      allCalls <- response.value.fieldAsObject("calls")
      callShards <- allCalls.fieldAsArray(jobKey.callFqn)
      onlyShardElement <- callShards.exactlyOneValueAsObject
      _ <- onlyShardElement.checkFieldValue("shardIndex", jobKey.index.getOrElse(-1).toString)
      callCaching <- onlyShardElement.fieldAsObject(CallMetadataKeys.CallCaching)
    } yield CallCachingMetadataJson(callCaching)
  }

  def extractHashes(callCachingMetadataJson: CallCachingMetadataJson): ErrorOr[Map[String, String]] = {
    def processField(keyPrefix: String)(fieldValue: (String, JsValue)): ErrorOr[Map[String, String]] = fieldValue match {
      case (key, hashString: JsString) => Map(keyPrefix + key -> hashString.value).validNel
      case (key, subObject: JsObject) => extractHashEntries(key + ":")(subObject)
      case other => s"Cannot extract hashes. Expected JsString or JsObject but got $other".invalidNel
    }

    def extractHashEntries(keyPrefix: String)(jsObject: JsObject): ErrorOr[Map[String, String]] = {
      val traversed = jsObject.fields.toList.traverse(processField(keyPrefix))
      traversed.map(_.flatten.toMap)
    }

    callCachingMetadataJson.value.fieldAsObject("hashes") flatMap extractHashEntries("")
  }

  def diffHashEvents(hashesA: Map[String, String], hashesB: Map[String, String]): Map[String, (Option[String], Option[String])] = {

    val hashesInANotMatchedInB: Map[String, (Option[String], Option[String])] = hashesA collect {
      case (key, value) if hashesB.get(key) != Option(value) => key -> (Option(value), hashesB.get(key))
    }

    val hashesUniqueToB: Map[String, String] = hashesB.filterNot({ case (k, _) => hashesA.keySet.contains(k) })

    hashesInANotMatchedInB ++ hashesUniqueToB.map {
      case (key, value) => key -> (None, Option(value))
    }
  }

  implicit class EnhancedJsObject(val jsObject: JsObject) extends AnyVal {
    def getField(field: String): ErrorOr[JsValue] = jsObject.fields.get(field).toErrorOr("No value provided")
    def fieldAsObject(field: String): ErrorOr[JsObject] = jsObject.getField(field) flatMap { _.mapToJsObject }
    def fieldAsArray(field: String): ErrorOr[JsArray] = jsObject.getField(field) flatMap { _.mapToJsArray }
    def checkFieldValue(field: String, expectation: String): ErrorOr[Unit] = jsObject.getField(field) flatMap {
      case v: JsValue if v.toString == expectation => ().validNel
      case other => s"Unexpected metadata field '$field'. Expected '$expectation' but got ${other.toString}".invalidNel
    }
  }

  implicit class EnhancedJsArray(val jsArray: JsArray) extends AnyVal {
    def exactlyOneValueAsObject: ErrorOr[JsObject] = if (jsArray.elements.size == 1) {
      jsArray.elements.head.mapToJsObject
    } else {
      s"Expected exactly 1 array element but got ${}".invalidNel
    }
  }

  implicit class EnhancedJsValue(val jsValue: JsValue) extends AnyVal {
    def mapToJsObject: ErrorOr[JsObject] = jsValue match {
      case obj: JsObject => obj.validNel
      case other => s"Invalid value type. Expected JsObject but got ${other.getClass.getSimpleName}: ${other.prettyPrint}".invalidNel
    }
    def mapToJsArray: ErrorOr[JsArray] = jsValue match {
      case arr: JsArray => arr.validNel
      case other => s"Invalid value type. Expected JsArray but got ${other.getClass.getSimpleName}: ${other.prettyPrint}".invalidNel
    }
  }
}
