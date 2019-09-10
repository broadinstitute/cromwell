package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.testkit.{ImplicitSender, TestFSMRef, TestProbe}
import cats.data.NonEmptyList
import cromwell.core._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffQueryParameter.CallCacheDiffQueryCall
import cromwell.services.metadata.MetadataService.GetMetadataAction
import cromwell.services.metadata.{MetadataService, _}
import cromwell.services.metadata.impl.builder.MetadataBuilderActor
import cromwell.services.metadata.impl.builder.MetadataBuilderActor.{BuiltMetadataResponse, FailedMetadataResponse}
import org.scalatest.concurrent.Eventually
import org.scalatest.{FlatSpecLike, Matchers}
import spray.json.JsObject

class CallCacheDiffActorSpec extends TestKitSuite with FlatSpecLike with Matchers with ImplicitSender with Eventually {

  behavior of "CallCacheDiffActor"

  val workflowIdA = WorkflowId.fromString("971652a6-139c-4ef3-96b5-aeb611a40dbf")
  val workflowIdB = WorkflowId.fromString("bb85b3ec-e179-4f12-b90f-5191216da598")

  val callFqnA = "callFqnA"
  val callFqnB = "callFqnB"

  val metadataJobKeyA = Option(MetadataJobKey(callFqnA, Option(1), 1))
  val metadataJobKeyB = Option(MetadataJobKey(callFqnB, None, 1))

  val callA = CallCacheDiffQueryCall(workflowIdA, callFqnA, Option(1))
  val callB = CallCacheDiffQueryCall(workflowIdB, callFqnB, None)

  val queryA = MetadataQuery(
    workflowId = workflowIdA,
    jobKey = Option(MetadataQueryJobKey(callFqnA, Option(1), None)),
    key = None,
    includeKeysOption = Option(NonEmptyList.of("callCaching", "executionStatus")),
    excludeKeysOption = Option(NonEmptyList.of("callCaching:hitFailures")),
    expandSubWorkflows = false
  )

  val queryB = MetadataQuery(
    workflowId = workflowIdB,
    jobKey = Option(MetadataQueryJobKey(callFqnB, None, None)),
    key = None,
    includeKeysOption = Option(NonEmptyList.of("callCaching", "executionStatus")),
    excludeKeysOption = Option(NonEmptyList.of("callCaching:hitFailures")),
    expandSubWorkflows = false
  )

  val eventsA = List(
      MetadataEvent(MetadataKey(workflowIdA, metadataJobKeyA, "executionStatus"), MetadataValue("Done")),
      MetadataEvent(MetadataKey(workflowIdA, metadataJobKeyA, "callCaching:allowResultReuse"), MetadataValue(true)),
      MetadataEvent(MetadataKey(workflowIdA, metadataJobKeyA, "callCaching:hashes:hash in only in A"), MetadataValue("hello from A")),
      MetadataEvent(MetadataKey(workflowIdA, metadataJobKeyA, "callCaching:hashes:hash in A and B with same value"), MetadataValue("we are thinking the same thought")),
      MetadataEvent(MetadataKey(workflowIdA, metadataJobKeyA, "callCaching:hashes:hash in A and B with different value"), MetadataValue("I'm the hash for A !"))
  )
  val workflowMetadataA: JsObject = MetadataBuilderActor.workflowMetadataResponse(workflowIdA, eventsA, includeCallsIfEmpty = false, Map.empty)
  val responseForA = BuiltMetadataResponse(MetadataService.GetMetadataAction(queryA), workflowMetadataA)

  val eventsB = List(
    MetadataEvent(MetadataKey(workflowIdB, metadataJobKeyB, "executionStatus"), MetadataValue("Failed")),
    MetadataEvent(MetadataKey(workflowIdB, metadataJobKeyB, "callCaching:allowResultReuse"), MetadataValue(false)),
    MetadataEvent(MetadataKey(workflowIdB, metadataJobKeyB, "callCaching:hashes:hash in only in B"), MetadataValue("hello from B")),
    MetadataEvent(MetadataKey(workflowIdB, metadataJobKeyB, "callCaching:hashes:hash in A and B with same value"), MetadataValue("we are thinking the same thought")),
    MetadataEvent(MetadataKey(workflowIdB, metadataJobKeyB, "callCaching:hashes:hash in A and B with different value"), MetadataValue("I'm the hash for B !"))
  )
  val workflowMetadataB: JsObject = MetadataBuilderActor.workflowMetadataResponse(workflowIdB, eventsB, includeCallsIfEmpty = false, Map.empty)
  val responseForB = BuiltMetadataResponse(MetadataService.GetMetadataAction(queryB), workflowMetadataB)

  it should "send correct queries to MetadataService when receiving a CallCacheDiffRequest" in {
    val mockServiceRegistryActor = TestProbe()
    val actor = TestFSMRef(new CallCacheDiffActor(mockServiceRegistryActor.ref))

    actor ! CallCacheDiffQueryParameter(callA, callB)

    mockServiceRegistryActor.expectMsg(GetMetadataAction(queryA))
    mockServiceRegistryActor.expectMsg(GetMetadataAction(queryB))

    system.stop(actor)
  }

  it should "save response for callA and wait for callB" in {
    val mockServiceRegistryActor = TestProbe()
    val actor = TestFSMRef(new CallCacheDiffActor(mockServiceRegistryActor.ref))

    actor.setState(WaitingForMetadata, CallCacheDiffWithRequest(queryA, queryB, None, None, self))

    actor ! responseForA

    eventually {
      actor.stateData shouldBe CallCacheDiffWithRequest(queryA, queryB, Some(WorkflowMetadataJson(workflowMetadataA)), None, self)
      actor.stateName shouldBe WaitingForMetadata
    }

    system.stop(actor)
  }

  it should "save response for callB and wait for callA" in {
    val mockServiceRegistryActor = TestProbe()
    val actor = TestFSMRef(new CallCacheDiffActor(mockServiceRegistryActor.ref))

    actor.setState(WaitingForMetadata, CallCacheDiffWithRequest(queryA, queryB, None, None, self))

    actor ! responseForB

    eventually {
      actor.stateData shouldBe CallCacheDiffWithRequest(queryA, queryB, None, Some(WorkflowMetadataJson(workflowMetadataB)), self)
      actor.stateName shouldBe WaitingForMetadata
    }

    system.stop(actor)
  }

  it should "build the response when receiving response for A and already has B" in {
    val mockServiceRegistryActor = TestProbe()
    val actor = TestFSMRef(new CallCacheDiffActor(mockServiceRegistryActor.ref))
    watch(actor)

    actor.setState(WaitingForMetadata, CallCacheDiffWithRequest(queryA, queryB, None, Some(WorkflowMetadataJson(workflowMetadataB)), self))

    actor ! responseForA

    expectMsgClass(classOf[CallCacheDiffActorResponse])
    expectTerminated(actor)
  }

  it should "build the response when receiving response for B and already has A" in {
    val mockServiceRegistryActor = TestProbe()
    val actor = TestFSMRef(new CallCacheDiffActor(mockServiceRegistryActor.ref))
    watch(actor)

    actor.setState(WaitingForMetadata, CallCacheDiffWithRequest(queryA, queryB, Some(WorkflowMetadataJson(workflowMetadataA)), None, self))

    actor ! responseForB

    expectMsgClass(classOf[CallCacheDiffActorResponse])
    expectTerminated(actor)
  }

  val correctCallCacheDiff = {
    import spray.json._

    s"""
       |{
       |   "callA":{
       |      "executionStatus": "Done",
       |      "allowResultReuse": true,
       |      "callFqn": "callFqnA",
       |      "jobIndex": 1,
       |      "workflowId": "971652a6-139c-4ef3-96b5-aeb611a40dbf"
       |   },
       |   "callB":{
       |      "executionStatus": "Failed",
       |      "allowResultReuse": false,
       |      "callFqn": "callFqnB",
       |      "jobIndex": -1,
       |      "workflowId": "bb85b3ec-e179-4f12-b90f-5191216da598"
       |   },
       |   "hashDifferential":[
       |      {
       |         "hashKey": "hash in only in A",
       |         "callA":"hello from A",
       |         "callB":null
       |      },
       |      {
       |         "hashKey": "hash in A and B with different value",
       |         "callA":"I'm the hash for A !",
       |         "callB":"I'm the hash for B !"
       |      },
       |      {
       |         "hashKey": "hash in only in B",
       |         "callA":null,
       |         "callB":"hello from B"
       |      }
       |   ]
       |}
       """.stripMargin.parseJson.asJsObject
  }

  it should "build a correct response" in {
    import spray.json._
    import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffActorJsonFormatting.successfulResponseJsonFormatter

    val mockServiceRegistryActor = TestProbe()
    val actor = TestFSMRef(new CallCacheDiffActor(mockServiceRegistryActor.ref))
    watch(actor)
    actor.setState(WaitingForMetadata, CallCacheDiffWithRequest(queryA, queryB, None, None, self))

    actor ! responseForB
    actor ! responseForA

    expectMsgPF() {
      case r: SuccessfulCallCacheDiffResponse =>
        withClue(s"""
             |Expected:
             |${correctCallCacheDiff.prettyPrint}
             |
             |Actual:
             |${r.toJson.prettyPrint}""".stripMargin) {
          r.toJson should be(correctCallCacheDiff)
        }
      case other => fail(s"Expected SuccessfulCallCacheDiffResponse but got $other")
    }
    expectTerminated(actor)
  }

  it should "build a correct response from multiple attempts" in {
    import spray.json._
    import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffActorJsonFormatting.successfulResponseJsonFormatter

    val mockServiceRegistryActor = TestProbe()
    val actor = TestFSMRef(new CallCacheDiffActor(mockServiceRegistryActor.ref))
    watch(actor)
    actor.setState(WaitingForMetadata, CallCacheDiffWithRequest(queryA, queryB, None, None, self))

    // Create a set of "failed" events for attempt 1:
    val eventsAAttempt1 = List(
      MetadataEvent(MetadataKey(workflowIdA, metadataJobKeyA, "executionStatus"), MetadataValue("Failed")),
      MetadataEvent(MetadataKey(workflowIdA, metadataJobKeyA, "callCaching:allowResultReuse"), MetadataValue(false)),
      MetadataEvent(MetadataKey(workflowIdA, metadataJobKeyA, "callCaching:hashes:hash in only in A"), MetadataValue("ouch!")),
      MetadataEvent(MetadataKey(workflowIdA, metadataJobKeyA, "callCaching:hashes:hash in A and B with same value"), MetadataValue("ouch!")),
      MetadataEvent(MetadataKey(workflowIdA, metadataJobKeyA, "callCaching:hashes:hash in A and B with different value"), MetadataValue("ouch!"))
    )
    // And update the old "eventsA" to represent attempt 2:
    val eventsAAttempt2 = eventsA.map(event => event.copy(key = event.key.copy(jobKey = event.key.jobKey.map(_.copy(attempt = 2)))))

    val modifiedEventsA = eventsAAttempt1 ++ eventsAAttempt2

    val workflowMetadataA: JsObject = MetadataBuilderActor.workflowMetadataResponse(workflowIdA, modifiedEventsA, includeCallsIfEmpty = false, Map.empty)
    val responseForA = BuiltMetadataResponse(MetadataService.GetMetadataAction(queryA), workflowMetadataA)


    actor ! responseForB
    actor ! responseForA

    expectMsgPF() {
      case r: SuccessfulCallCacheDiffResponse =>
        withClue(s"""
                    |Expected:
                    |${correctCallCacheDiff.prettyPrint}
                    |
                    |Actual:
                    |${r.toJson.prettyPrint}""".stripMargin) {
          r.toJson should be(correctCallCacheDiff)
        }
      case other =>
        expectTerminated(actor)
        fail(s"Expected SuccessfulCallCacheDiffResponse but got $other")
    }

    expectTerminated(actor)
  }

  it should "fail properly" in {
    import scala.concurrent.duration._
    import scala.language.postfixOps

    val mockServiceRegistryActor = TestProbe()
    val actor = TestFSMRef(new CallCacheDiffActor(mockServiceRegistryActor.ref))
    watch(actor)
    val exception = new Exception("Query lookup failed - but it's ok ! this is a test !")
    val responseA = FailedMetadataResponse(GetMetadataAction(queryA), exception)

    actor.setState(WaitingForMetadata, CallCacheDiffWithRequest(queryA, queryB, None, None, self))

    actor ! responseA

    expectMsgPF(1 second) {
      case FailedCallCacheDiffResponse(e: Throwable) =>
        e.getMessage shouldBe "Query lookup failed - but it's ok ! this is a test !"
    }

    expectTerminated(actor)
  }

  it should "respond with an appropriate error if calls' hashes are missing" in {
    testExpectedErrorForModifiedMetadata(
      metadataFilter = _.key.key.contains("hashes"),
      error = s"""Failed to calculate diff for call A and call B:
                 |Failed to extract relevant metadata for call A (971652a6-139c-4ef3-96b5-aeb611a40dbf / callFqnA:1) (reason 1 of 1): No 'hashes' field found
                 |Failed to extract relevant metadata for call B (bb85b3ec-e179-4f12-b90f-5191216da598 / callFqnB:-1) (reason 1 of 1): No 'hashes' field found""".stripMargin
    )
  }

  it should "respond with an appropriate error if both calls are missing" in {
    testExpectedErrorForModifiedMetadata(
      metadataFilter = _.key.jobKey.nonEmpty,
      error = s"""Failed to calculate diff for call A and call B:
                 |Failed to extract relevant metadata for call A (971652a6-139c-4ef3-96b5-aeb611a40dbf / callFqnA:1) (reason 1 of 1): No 'calls' field found
                 |Failed to extract relevant metadata for call B (bb85b3ec-e179-4f12-b90f-5191216da598 / callFqnB:-1) (reason 1 of 1): No 'calls' field found""".stripMargin
    )
  }

  def testExpectedErrorForModifiedMetadata(metadataFilter: MetadataEvent => Boolean, error: String) = {
    import scala.concurrent.duration._
    import scala.language.postfixOps

    val mockServiceRegistryActor = TestProbe()
    val actor = TestFSMRef(new CallCacheDiffActor(mockServiceRegistryActor.ref))
    watch(actor)

    def getModifiedResponse(workflowId: WorkflowId, query: MetadataQuery, events: Seq[MetadataEvent]): BuiltMetadataResponse = {
      val modifiedEvents = events.filterNot(metadataFilter) // filters out any "call" level metadata
      val modifiedWorkflowMetadata = MetadataBuilderActor.workflowMetadataResponse(workflowId, modifiedEvents, includeCallsIfEmpty = false, Map.empty)
      BuiltMetadataResponse(MetadataService.GetMetadataAction(query), modifiedWorkflowMetadata)
    }

    actor.setState(WaitingForMetadata, CallCacheDiffWithRequest(queryA, queryB, None, None, self))

    actor ! getModifiedResponse(workflowIdA, queryA, eventsA)
    actor ! getModifiedResponse(workflowIdB, queryB, eventsB)

    expectMsgPF(1 second) {
      case FailedCallCacheDiffResponse(e) => e.getMessage shouldBe error
    }
    expectTerminated(actor)
  }

}
