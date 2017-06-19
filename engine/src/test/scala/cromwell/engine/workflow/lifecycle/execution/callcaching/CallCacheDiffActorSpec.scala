package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.testkit.{ImplicitSender, TestFSMRef, TestProbe}
import cats.data.NonEmptyList
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffActor.{CallCacheDiffWithRequest, WaitingForMetadata}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffQueryParameter.CallCacheDiffQueryCall
import cromwell.services.metadata.MetadataService.{GetMetadataQueryAction, MetadataLookupResponse, MetadataServiceKeyLookupFailed}
import cromwell.services.metadata._
import cromwell.webservice.FailureResponse
import cromwell.webservice.PerRequest.RequestComplete
import org.scalatest.concurrent.Eventually
import org.scalatest.{FlatSpecLike, Matchers}
import spray.http.StatusCodes

class CallCacheDiffActorSpec extends TestKitSuite with FlatSpecLike with Matchers with ImplicitSender with Eventually {

  behavior of "CallCacheDiffActor"

  val workflowIdA = WorkflowId.fromString("971652a6-139c-4ef3-96b5-aeb611a40dbf")
  val workflowIdB = WorkflowId.fromString("bb85b3ec-e179-4f12-b90f-5191216da598")

  val callFqnA = "callFqnA"
  val callFqnB = "callFqnB"
  
  val metadataJobKeyA = Option(MetadataJobKey(callFqnA, Option(1), 1))
  val metadataJobKeyB = Option(MetadataJobKey(callFqnB, None, 1))
  
  val callA = CallCacheDiffQueryCall(workflowIdA.toString, callFqnA, Option(1))
  val callB = CallCacheDiffQueryCall(workflowIdB.toString, callFqnB, None)

  val queryA = MetadataQuery(
    workflowIdA,
    Option(MetadataQueryJobKey(callFqnA, Option(1), None)),
    None,
    Option(NonEmptyList.of("callCaching", "executionStatus")),
    None,
    expandSubWorkflows = false
  )

  val queryB = MetadataQuery(
    workflowIdB,
    Option(MetadataQueryJobKey(callFqnB, None, None)),
    None,
    Option(NonEmptyList.of("callCaching", "executionStatus")),
    None,
    expandSubWorkflows = false
  )
  
  val eventsA = List(
      MetadataEvent(MetadataKey(workflowIdA, metadataJobKeyA, "executionStatus"), MetadataValue("Done")),
      MetadataEvent(MetadataKey(workflowIdA, metadataJobKeyA, "callCaching:allowResultReuse"), MetadataValue(true)),
      MetadataEvent(MetadataKey(workflowIdA, metadataJobKeyA, "callCaching:hashes: hash in only in A"), MetadataValue("hello")),
      MetadataEvent(MetadataKey(workflowIdA, metadataJobKeyA, "callCaching:hashes: hash in A and B with same value"), MetadataValue(1)),
      MetadataEvent(MetadataKey(workflowIdA, metadataJobKeyA, "callCaching:hashes: hash in A and B with different value"), MetadataValue("I'm the hash for A !"))
  )

  val eventsB = List(
    MetadataEvent(MetadataKey(workflowIdB, metadataJobKeyB, "executionStatus"), MetadataValue("Failed")),
    MetadataEvent(MetadataKey(workflowIdB, metadataJobKeyB, "callCaching:allowResultReuse"), MetadataValue(false)),
    MetadataEvent(MetadataKey(workflowIdB, metadataJobKeyB, "callCaching:hashes: hash in only in B"), MetadataValue("hello")),
    MetadataEvent(MetadataKey(workflowIdB, metadataJobKeyB, "callCaching:hashes: hash in A and B with same value"), MetadataValue(1)),
    MetadataEvent(MetadataKey(workflowIdA, metadataJobKeyA, "callCaching:hashes: hash in A and B with different value"), MetadataValue("I'm the hash for B !"))
  )
  
  it should "send correct queries to MetadataService when receiving a CallCacheDiffRequest" in {
    val mockServiceRegistryActor = TestProbe()
    val actor = TestFSMRef(new CallCacheDiffActor(mockServiceRegistryActor.ref))

    actor ! CallCacheDiffQueryParameter(callA, callB)

    mockServiceRegistryActor.expectMsg(GetMetadataQueryAction(queryA))
    mockServiceRegistryActor.expectMsg(GetMetadataQueryAction(queryB))
  }

  it should "save response for callA and wait for callB" in {
    val mockServiceRegistryActor = TestProbe()
    val actor = TestFSMRef(new CallCacheDiffActor(mockServiceRegistryActor.ref))
    
    actor.setState(WaitingForMetadata, CallCacheDiffWithRequest(queryA, queryB, None, None, self))

    val response = MetadataLookupResponse(queryA, eventsA)
    actor ! response
   
    eventually {
      actor.stateData shouldBe CallCacheDiffWithRequest(queryA, queryB, Some(response), None, self)
      actor.stateName shouldBe WaitingForMetadata
    }
  }

  it should "save response for callB and wait for callA" in {
    val mockServiceRegistryActor = TestProbe()
    val actor = TestFSMRef(new CallCacheDiffActor(mockServiceRegistryActor.ref))

    actor.setState(WaitingForMetadata, CallCacheDiffWithRequest(queryA, queryB, None, None, self))

    val response = MetadataLookupResponse(queryB, eventsB)
    actor ! response

    eventually {
      actor.stateData shouldBe CallCacheDiffWithRequest(queryA, queryB, None, Some(response), self)
      actor.stateName shouldBe WaitingForMetadata
    }
  }

  it should "build the response when receiving response for A and already has B" in {
    val mockServiceRegistryActor = TestProbe()
    val actor = TestFSMRef(new CallCacheDiffActor(mockServiceRegistryActor.ref))
    watch(actor)
    val responseB = MetadataLookupResponse(queryB, eventsB)

    actor.setState(WaitingForMetadata, CallCacheDiffWithRequest(queryA, queryB, None, Option(responseB), self))

    actor ! MetadataLookupResponse(queryA, eventsA)

    expectMsgClass(classOf[RequestComplete[_]])
    expectTerminated(actor)
  }

  it should "build the response when receiving response for B and already has A" in {
    val mockServiceRegistryActor = TestProbe()
    val actor = TestFSMRef(new CallCacheDiffActor(mockServiceRegistryActor.ref))
    watch(actor)
    val responseA = MetadataLookupResponse(queryA, eventsA)

    actor.setState(WaitingForMetadata, CallCacheDiffWithRequest(queryA, queryB, Option(responseA), None, self))

    actor ! MetadataLookupResponse(queryB, eventsB)

    expectMsgClass(classOf[RequestComplete[_]])
    expectTerminated(actor)
  }

  it should "build a correct response" in {
    import cromwell.services.metadata.MetadataService.MetadataLookupResponse
    import cromwell.webservice.PerRequest.RequestComplete
    import cromwell.webservice.WorkflowJsonSupport._
    import spray.http.StatusCodes
    import spray.httpx.SprayJsonSupport._
    import spray.json._
    
    val mockServiceRegistryActor = TestProbe()
    val actor = TestFSMRef(new CallCacheDiffActor(mockServiceRegistryActor.ref))
    watch(actor)
    actor.setState(WaitingForMetadata, CallCacheDiffWithRequest(queryA, queryB, None, None, self))

    actor ! MetadataLookupResponse(queryB, eventsB)
    actor ! MetadataLookupResponse(queryA, eventsA)

    val expectedJson: JsObject =
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
         |         "hash in only in A":{  
         |            "callA":"hello",
         |            "callB":null
         |         }
         |      },
         |      {  
         |         "hash in A and B with different value":{  
         |            "callA":"I'm the hash for A !",
         |            "callB":"I'm the hash for B !"
         |         }
         |      },
         |      {  
         |         "hash in only in B":{  
         |            "callA":null,
         |            "callB":"hello"
         |         }
         |      }
         |   ]
         |}
       """.stripMargin.parseJson.asJsObject
    
    val expectedResponse = RequestComplete((StatusCodes.OK, expectedJson))
    
    expectMsg(expectedResponse)
    expectTerminated(actor)
  }

  it should "fail properly" in {
    import scala.concurrent.duration._
    import scala.language.postfixOps

    val mockServiceRegistryActor = TestProbe()
    val actor = TestFSMRef(new CallCacheDiffActor(mockServiceRegistryActor.ref))
    watch(actor)
    val exception = new Exception("Query lookup failed - but it's ok ! this is a test !")
    val responseA = MetadataServiceKeyLookupFailed(queryA, exception)

    actor.setState(WaitingForMetadata, CallCacheDiffWithRequest(queryA, queryB, None, None, self))

    actor ! responseA

    expectMsgPF(1 second) {
      case RequestComplete((StatusCodes.InternalServerError, response: FailureResponse)) =>
        response.status shouldBe "error"
        response.message shouldBe "Query lookup failed - but it's ok ! this is a test !"
    }
    
    expectTerminated(actor)
  }

}
