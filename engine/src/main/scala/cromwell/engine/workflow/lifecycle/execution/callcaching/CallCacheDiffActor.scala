package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{ActorRef, LoggingFSM, Props}
import cats.data.NonEmptyList
import cats.instances.list._
import cats.syntax.apply._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.exception.AggregatedMessageException
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffActor.{CallCacheDiffActorData, _}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffQueryParameter.CallCacheDiffQueryCall
import cromwell.services.metadata.MetadataService.GetMetadataAction
import cromwell.services.metadata._
import cromwell.services.metadata.impl.builder.MetadataBuilderActor.{BuiltMetadataResponse, FailedMetadataResponse}
import spray.json.{JsArray, JsBoolean, JsNumber, JsObject, JsString, JsValue}

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
    case Event(FailedMetadataResponse(_, failure), data: CallCacheDiffWithRequest) =>
      data.replyTo ! FailedCallCacheDiffResponse(failure)
      context stop self
      stay()
  }

  whenUnhandled {
    case Event(oops, oopsData) =>
      log.error(s"Programmer Error: Unexpected event received by ${this.getClass.getSimpleName}: $oops / $oopsData (in state $stateName)")
      stay()

  }

  private def buildDiffAndRespond(queryA: MetadataQuery,
                                  queryB: MetadataQuery,
                                  responseA: WorkflowMetadataJson,
                                  responseB: WorkflowMetadataJson,
                                  replyTo: ActorRef) = {

    def describeCallFromQuery(query: MetadataQuery): String = s"${query.workflowId} / ${query.jobKey.map(_.callFqn).getOrElse("<<CallMissing>>")}:${query.jobKey.map(_.index.getOrElse(-1)).getOrElse("<<CallMissing>>")}"

    val callACachingMetadata = extractCallMetadata(queryA, responseA).contextualizeErrors(s"extract relevant metadata for call A (${describeCallFromQuery(queryA)})")
    val callBCachingMetadata = extractCallMetadata(queryB, responseB).contextualizeErrors(s"extract relevant metadata for call B (${describeCallFromQuery(queryB)})")

    val response = (callACachingMetadata, callBCachingMetadata) flatMapN { case (callA, callB) =>

        val callADetails = extractCallDetails(queryA, callA)
        val callBDetails = extractCallDetails(queryB, callB)

      (callADetails, callBDetails) mapN { (cad, cbd) =>
        val callAHashes = callA.callCachingMetadataJson.hashes
        val callBHashes = callB.callCachingMetadataJson.hashes

        SuccessfulCallCacheDiffResponse(cad, cbd, calculateHashDifferential(callAHashes, callBHashes))
      }
    } valueOr {
      e => FailedCallCacheDiffResponse(AggregatedMessageException("Failed to calculate diff for call A and call B", e.toList))
    }

    replyTo ! response

    context stop self
    stay()
  }
}


object CallCacheDiffActor {

  final case class CachedCallNotFoundException(message: String) extends Exception {
    override def getMessage = message
  }

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

  case class FailedCallCacheDiffResponse(reason: Throwable) extends CallCacheDiffActorResponse
  final case class SuccessfulCallCacheDiffResponse(callA: CallDetails, callB: CallDetails, hashDifferential: List[HashDifference]) extends CallCacheDiffActorResponse
  def props(serviceRegistryActor: ActorRef) = Props(new CallCacheDiffActor(serviceRegistryActor)).withDispatcher(EngineDispatcher)

  final case class CallDetails(executionStatus: String, allowResultReuse: Boolean, callFqn: String, jobIndex: Int, workflowId: String)
  final case class HashDifference(hashKey: String, callA: Option[String], callB: Option[String])


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
  final case class CallMetadataJson(rawValue: JsObject, jobKey: MetadataQueryJobKey, callCachingMetadataJson: CallCachingMetadataJson)
  final case class CallCachingMetadataJson(rawValue: JsObject, hashes: Map[String, String])


  /*
   * Takes in the JsObject returned from a metadata query and filters out only the appropriate call's callCaching section
   */
  def extractCallMetadata(query: MetadataQuery, response: WorkflowMetadataJson): ErrorOr[CallMetadataJson] = {

    for {
      // Sanity Checks:
      _ <- response.value.checkFieldValue("id", s""""${query.workflowId}"""")
      jobKey <- query.jobKey.toErrorOr("Call is required in call cache diff query")

      // Unpack the JSON:
      allCalls <- response.value.fieldAsObject("calls")
      callShards <- allCalls.fieldAsArray(jobKey.callFqn)
      onlyShardElement <- callShards.elementWithHighestAttemptField
      _ <- onlyShardElement.checkFieldValue("shardIndex", jobKey.index.getOrElse(-1).toString)
      callCachingElement <- onlyShardElement.fieldAsObject(CallMetadataKeys.CallCaching)
      hashes <- extractHashes(callCachingElement)
    } yield CallMetadataJson(onlyShardElement, jobKey, CallCachingMetadataJson(callCachingElement, hashes))
  }

  def extractHashes(callCachingMetadataJson: JsObject): ErrorOr[Map[String, String]] = {
    def processField(keyPrefix: String)(fieldValue: (String, JsValue)): ErrorOr[Map[String, String]] = fieldValue match {
      case (key, hashString: JsString) => Map(keyPrefix + key -> hashString.value).validNel
      case (key, subObject: JsObject) => extractHashEntries(key + ":", subObject)
      case (key, otherValue) => s"Cannot extract hashes for $key. Expected JsString or JsObject but got ${otherValue.getClass.getSimpleName} $otherValue".invalidNel
    }

    def extractHashEntries(keyPrefix: String, jsObject: JsObject): ErrorOr[Map[String, String]] = {
      val traversed = jsObject.fields.toList.traverse(processField(keyPrefix))
      traversed.map(_.flatten.toMap)
    }

    for {
      hashesSection <- callCachingMetadataJson.fieldAsObject("hashes")
      entries <- extractHashEntries("", hashesSection)
    } yield entries
  }

  def calculateHashDifferential(hashesA: Map[String, String], hashesB: Map[String, String]): List[HashDifference] = {
    val hashesInANotMatchedInB: List[HashDifference] = hashesA.toList collect {
      case (key, value) if hashesB.get(key) != Option(value) => HashDifference(key, Option(value), hashesB.get(key))
    }
    val hashesUniqueToB: List[HashDifference] = hashesB.toList.collect {
      case (key, value) if !hashesA.keySet.contains(key) => HashDifference(key, None, Option(value))
    }
    hashesInANotMatchedInB ++ hashesUniqueToB
  }

  def extractCallDetails(query: MetadataQuery, callMetadataJson: CallMetadataJson): ErrorOr[CallDetails] = {
    val executionStatus = callMetadataJson.rawValue.fieldAsString("executionStatus")
    val allowResultReuse = callMetadataJson.callCachingMetadataJson.rawValue.fieldAsBoolean("allowResultReuse")

    (executionStatus, allowResultReuse) mapN { (es, arr) =>
      CallDetails(
        executionStatus = es.value,
        allowResultReuse = arr.value,
        callFqn = callMetadataJson.jobKey.callFqn,
        jobIndex = callMetadataJson.jobKey.index.getOrElse(-1),
        workflowId = query.workflowId.toString
      )
    }
  }

  implicit class EnhancedJsObject(val jsObject: JsObject) extends AnyVal {
    def getField(field: String): ErrorOr[JsValue] = jsObject.fields.get(field).toErrorOr(s"No '$field' field found")
    def fieldAsObject(field: String): ErrorOr[JsObject] = jsObject.getField(field) flatMap { _.mapToJsObject }
    def fieldAsArray(field: String): ErrorOr[JsArray] = jsObject.getField(field) flatMap { _.mapToJsArray }
    def fieldAsString(field: String): ErrorOr[JsString] = jsObject.getField(field) flatMap { _.mapToJsString }
    def fieldAsNumber(field: String): ErrorOr[JsNumber] = jsObject.getField(field) flatMap { _.mapToJsNumber }
    def fieldAsBoolean(field: String): ErrorOr[JsBoolean] = jsObject.getField(field) flatMap { _.mapToJsBoolean }
    def checkFieldValue(field: String, expectation: String): ErrorOr[Unit] = jsObject.getField(field) flatMap {
      case v: JsValue if v.toString == expectation => ().validNel
      case other => s"Unexpected metadata field '$field'. Expected '$expectation' but got ${other.toString}".invalidNel
    }
  }

  implicit class EnhancedJsArray(val jsArray: JsArray) extends AnyVal {

    def elementWithHighestAttemptField: ErrorOr[JsObject] = {
      def extractAttemptAndObject(value: JsValue): ErrorOr[(Int, JsObject)] = for {
        asObject <- value.mapToJsObject
        attempt <- asObject.fieldAsNumber("attempt")
      } yield (attempt.value.intValue(), asObject)

      def foldFunction(accumulator: ErrorOr[(Int, JsObject)], nextElement: JsValue): ErrorOr[(Int, JsObject)] = {
        (accumulator, extractAttemptAndObject(nextElement)) mapN { case ((previousHighestAttempt, previousJsObject), (nextAttempt, nextJsObject)) =>
          if (previousHighestAttempt > nextAttempt) {
            (previousHighestAttempt, previousJsObject)
          } else {
            (nextAttempt, nextJsObject)
          }
        }
      }

      for {
        attemptListNel <- NonEmptyList.fromList(jsArray.elements.toList).toErrorOr("Expected at least one attempt but found 0")
        highestAttempt <- attemptListNel.toList.foldLeft(extractAttemptAndObject(attemptListNel.head))(foldFunction)
      } yield highestAttempt._2
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
    def mapToJsString: ErrorOr[JsString] = jsValue match {
      case str: JsString => str.validNel
      case other => s"Invalid value type. Expected JsString but got ${other.getClass.getSimpleName}: ${other.prettyPrint}".invalidNel
    }
    def mapToJsBoolean: ErrorOr[JsBoolean] = jsValue match {
      case boo: JsBoolean => boo.validNel
      case other => s"Invalid value type. Expected JsBoolean but got ${other.getClass.getSimpleName}: ${other.prettyPrint}".invalidNel
    }
    def mapToJsNumber: ErrorOr[JsNumber] = jsValue match {
      case boo: JsNumber => boo.validNel
      case other => s"Invalid value type. Expected JsNumber but got ${other.getClass.getSimpleName}: ${other.prettyPrint}".invalidNel
    }
  }
}
