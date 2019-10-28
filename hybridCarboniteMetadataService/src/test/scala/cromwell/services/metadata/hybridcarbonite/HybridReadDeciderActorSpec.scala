package cromwell.services.metadata.hybridcarbonite

import akka.actor.ActorSystem
import akka.testkit.{TestFSMRef, TestProbe}
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.services.SuccessfulMetadataJsonResponse
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataArchiveStatus._
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata.hybridcarbonite.HybridReadDeciderActor._
import org.scalatest.{FlatSpecLike, Matchers}
import spray.json._

class HybridReadDeciderActorSpec extends TestKitSuite("HybridReadDeciderActorSpec") with FlatSpecLike with Matchers {

  behavior of "HybridReadDeciderActor"

  implicit val as: ActorSystem = system

  it should "pass messages correctly for not-yet-archived workflow queries" in { routingTest(Unarchived, shouldChooseCarboniter = false) }
  it should "pass messages correctly for archived workflow queries" in { routingTest(Archived, shouldChooseCarboniter = true) }
  it should "pass messages correctly for failed-archived workflow queries" in { routingTest(ArchiveFailed, shouldChooseCarboniter = false) }

  def routingTest(archiveStatusToReturn: MetadataArchiveStatus, shouldChooseCarboniter: Boolean) = {
    val sampleWorkflowId = WorkflowId.randomId()

    val client = TestProbe("client")
    val classicMetadataActor = TestProbe("classic")
    val carboniteMetadataActor = TestProbe("carboniter")

    val decidedUponActor = if (shouldChooseCarboniter) carboniteMetadataActor else classicMetadataActor

    val hrda = TestFSMRef.apply(new HybridReadDeciderActor(classicMetadataActor.ref, carboniteMetadataActor.ref))
    watch(hrda)

    hrda.stateName should be(Pending)

    val incomingQuery = WorkflowOutputs(sampleWorkflowId)

    // The HybridReadDeciderActor gets a workflow-specific query and looks up its archive status:
    client.send(hrda, incomingQuery)
    classicMetadataActor.expectMsg(QueryForWorkflowsMatchingParameters(Vector("Id" -> sampleWorkflowId.toString)))
    hrda.stateName should be(RequestingMetadataArchiveStatus)

    // The classic metadata actor lets us know that the workflow is unarchived
    classicMetadataActor.send(hrda, WorkflowQuerySuccess(WorkflowQueryResponse(Seq(WorkflowQueryResult(null, null, null, null, null, null, null, null, null, archiveStatusToReturn)), 1), None))

    // We expect the classic metadata service to receive the original query now:
    decidedUponActor.expectMsg(incomingQuery)
    hrda.stateName should be(WaitingForMetadataResponse)

    // Send a basic response:
    val responseJsonString =
      s"""{
        |  "outputs": {},
        |  "id": "${sampleWorkflowId.toString}"
        |}""".stripMargin
    val responseJson = responseJsonString.parseJson.asJsObject
    val response = SuccessfulMetadataJsonResponse(incomingQuery, responseJson)
    decidedUponActor.send(hrda, response)

    client.expectMsg(response)

    // Expect the HRDA to shut itself down now that its job is done:
    expectTerminated(hrda)
    client.msgAvailable should be(false)
    classicMetadataActor.msgAvailable should be(false)
    carboniteMetadataActor.msgAvailable should be(false)
  }

  it should "go straight to the classic metadata service for summary table searches" in {

    val client = TestProbe("client")
    val classicMetadataActor = TestProbe("classic")
    val carboniteMetadataActor = TestProbe("carboniter")

    val hrda = TestFSMRef.apply(new HybridReadDeciderActor(classicMetadataActor.ref, carboniteMetadataActor.ref))
    watch(hrda)

    hrda.stateName should be(Pending)

    val queryMsg = QueryForWorkflowsMatchingParameters(Vector("Includekey" -> "blah"))
    client.send(hrda, queryMsg)
    classicMetadataActor.expectMsg(queryMsg)
    hrda.stateName should be(WaitingForMetadataResponse)

    val response = WorkflowQuerySuccess(WorkflowQueryResponse(Seq(
      WorkflowQueryResult("id1", null, null, null, null, null, null, null, null, null),
      WorkflowQueryResult("id2", null, null, null, null, null, null, null, null, null)
    ), 2), None)
    classicMetadataActor.send(hrda, response)
    client.expectMsg(response)

    // Expect the HRDA to shut itself down now that its job is done:
    expectTerminated(hrda)
    client.msgAvailable should be(false)
    classicMetadataActor.msgAvailable should be(false)
    carboniteMetadataActor.msgAvailable should be(false)
  }

}
