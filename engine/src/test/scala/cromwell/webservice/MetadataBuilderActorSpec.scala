package cromwell.webservice

import java.sql.Timestamp

import akka.actor.ActorSystem
import akka.testkit._
import cromwell.core.{KnowsWhatTimeItIs, WorkflowId}
import cromwell.services.MetadataServiceActor._
import cromwell.services._
import cromwell.webservice.PerRequest.RequestComplete
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{BeforeAndAfterAll, FlatSpecLike, Matchers}
import spray.http.{StatusCode, StatusCodes}
import spray.json._

import scala.concurrent.duration._
import scala.language.postfixOps

class MetadataBuilderActorSpec extends TestKit(ActorSystem("Metadata"))
  with FlatSpecLike with Matchers with TableDrivenPropertyChecks with ImplicitSender with BeforeAndAfterAll with KnowsWhatTimeItIs {

  behavior of "MetadataParser"

  val defaultTimeout = 100 millis
  val mockServiceRegistry = TestProbe()
  val parentProbe = TestProbe()

  def assertMetadataResponse(action: MetadataServiceAction,
                             queryReply: MetadataQuery,
                             events: Seq[MetadataEvent],
                             expectedRes: String) = {
    val metadataBuilder = TestActorRef(MetadataBuilderActor.props(mockServiceRegistry.ref), parentProbe.ref, "MetadataActor")
    metadataBuilder ! action // Ask for everything
    mockServiceRegistry.expectMsg(defaultTimeout, action) // TestActor runs on CallingThreadDispatcher
    mockServiceRegistry.reply(MetadataLookupResponse(queryReply, events))

    parentProbe.expectMsgPF(defaultTimeout) {
      case response: RequestComplete[(StatusCode, JsObject)] @unchecked =>
        response.response._1 shouldBe StatusCodes.OK
        response.response._2 shouldBe expectedRes.parseJson
    }
  }

  it should "build workflow scope tree from metadata events" in {
    def makeEvent(workflow: WorkflowId, key: Option[MetadataJobKey]) = {
      MetadataEvent(MetadataKey(workflow, key, "NOT_CHECKED"), MetadataValue("NOT_CHECKED"), currentTime)
    }

    val workflowA = WorkflowId.randomId()

    val workflowACalls = List(
      Option(new MetadataJobKey("callB", Some(1), 3)),
      Option(new MetadataJobKey("callB", None, 1)),
      Option(new MetadataJobKey("callB", Some(1), 2)),
      Option(new MetadataJobKey("callA", None, 1)),
      Option(new MetadataJobKey("callB", Some(1), 1)),
      Option(new MetadataJobKey("callB", Some(0), 1)),
      None
    )
    val workflowAEvents = workflowACalls map { makeEvent(workflowA, _) }

    val expectedRes =
      s"""{
        |  "$workflowA": {
        |    "calls": {
        |      "callB": [{
        |        "attempt": 1,
        |        "NOT_CHECKED": "NOT_CHECKED",
        |        "shardIndex": -1
        |      }, {
        |        "attempt": 1,
        |        "NOT_CHECKED": "NOT_CHECKED",
        |        "shardIndex": 0
        |      }, {
        |        "attempt": 1,
        |        "NOT_CHECKED": "NOT_CHECKED",
        |        "shardIndex": 1
        |      }, {
        |        "attempt": 2,
        |        "NOT_CHECKED": "NOT_CHECKED",
        |        "shardIndex": 1
        |      }, {
        |        "attempt": 3,
        |        "NOT_CHECKED": "NOT_CHECKED",
        |        "shardIndex": 1
        |      }],
        |      "callA": [{
        |        "attempt": 1,
        |        "NOT_CHECKED": "NOT_CHECKED",
        |        "shardIndex": -1
        |      }]
        |    },
        |    "NOT_CHECKED": "NOT_CHECKED"
        |  }
        |}""".stripMargin

    val mdQuery = MetadataQuery(workflowA, None, None)
    val queryAction = GetMetadataQueryAction(mdQuery)
    assertMetadataResponse(queryAction, mdQuery, workflowAEvents, expectedRes)
  }

  type EventBuilder = (String, String, Timestamp)

  def assertMetadataKeyStructure(eventList: List[EventBuilder], expectedJson: String) = {
    val workflow = WorkflowId.randomId()

    def makeEvent(workflow: WorkflowId)(key: String, value: String, timestamp: Timestamp) = {
      MetadataEvent(MetadataKey(workflow, None, key), MetadataValue(value), timestamp)
    }

    val events = eventList map Function.tupled(makeEvent(workflow))
    val expectedRes = s"""{"${workflow.id.toString}": { "calls": {}, $expectedJson } }"""

    val mdQuery = MetadataQuery(workflow, None, None)
    val queryAction = GetAllMetadataAction(workflow)
    assertMetadataResponse(queryAction, mdQuery, events, expectedRes)
  }

  it should "keep event with later timestamp for the same key in metadata" in {
    val eventBuilderList = List(
      ("a", "aLater", Timestamp.valueOf("2000-01-02 12:00:00")),
      ("a", "a", Timestamp.valueOf("2000-01-01 12:00:00"))
    )
    val expectedRes =
      """"a": "aLater"""".stripMargin

    assertMetadataKeyStructure(eventBuilderList, expectedRes)
  }

  it should "build JSON object structure from dotted key syntax" in {
    val eventBuilderList = List(
      ("a:b:c", "abc", currentTime),
      ("b:a", "ba", currentTime),
      ("c", "c", currentTime)
    )

    val expectedRes =
      """"a": {
        |    "b": {
        |      "c": "abc"
        |    }
        |  },
        |  "b": {
        |    "a": "ba"
        |  },
        |  "c": "c"""".stripMargin

    assertMetadataKeyStructure(eventBuilderList, expectedRes)
  }


  it should "build lexigraphically sorted JSON list structure from dotted key syntax" in {
    val eventBuilderList = List(
      ("l1[2]", "l12", currentTime),
      ("l1[k]", "l1k", currentTime),
      ("l1[3]", "l13", currentTime),
      ("l1[i]", "l1i", currentTime),
      ("l1[10]", "l110", currentTime),
      ("l1[j]", "l1j", currentTime)
    )

    val expectedRes =
      """"l1": [
        |    "l110", "l12", "l13", "l1i", "l1j", "l1k"
        |  ]""".stripMargin

    assertMetadataKeyStructure(eventBuilderList, expectedRes)
  }

  it should "override elements with same index in a list if they can't be merged together" in {
    val eventBuilderList = List(
      ("l1[1]", "a", currentTime),
      ("l1[2]", "a", currentTime),
      ("l1[3]", "a", currentTime),
      ("l1[2]", "b", currentTime),
      ("l1[3]", "b", currentTime),
      ("l1[3]", "c", currentTime)
    )

    val expectedRes =
      """"l1": [
        |    "a", "b", "c"
        |  ]""".stripMargin

    assertMetadataKeyStructure(eventBuilderList, expectedRes)
  }

  it should "nest lists and objects together and respect ordering" in {
    val eventBuilderList = List(
      ("l1[0]:l11[0]:l111[2]", "l10l110l1112", currentTime),
      ("l1[0]:l11[1]:l112[1]", "l10l110l1121", currentTime),
      ("l1[0]:l11[1]:l112[0]", "l10l110l1120", currentTime),
      ("l1[1]:a:l12[2]:b", "l11al122b", currentTime),
      ("l1[0]:l11[0]:l111[0]", "l10l110l1110", currentTime),
      ("l1[0]:l11[0]:l111[1]", "l10l110l1111", currentTime),
      ("l1[0]:l11[2]:l121[0]:a:b", "l10l111l1210ab", currentTime),
      ("l1[1]:a:l12[0]:b", "l11al120b", currentTime),
      ("l1[0]:l11[1]:l112[2]", "l10l110l1122", currentTime),
      ("l1[1]:a:l12[1]", "l11al121", currentTime),
      ("l1[0]:l11[2]:l121[0]:a:c", "l10l111l1210ac", currentTime)
    )

    val expectedRes =
      """"l1": [
        |      { "l11": [
        |      { "l111": ["l10l110l1110", "l10l110l1111", "l10l110l1112"] },
        |      { "l112": ["l10l110l1120", "l10l110l1121", "l10l110l1122"] },
        |      { "l121": [
        |        {
        |          "a": {
        |            "c": "l10l111l1210ac",
        |            "b": "l10l111l1210ab"
        |          }
        |        }]
        |      }]
        |    }, {
        |      "a": {
        |        "l12": [{
        |          "b": "l11al120b"
        |        }, "l11al121", {
        |          "b": "l11al122b"
        |        }]
        |      }
        |    }]""".stripMargin

    assertMetadataKeyStructure(eventBuilderList, expectedRes)
  }

  it should "override json values if they can't be merged" in {
    val kv = ("key", "value", currentTime)
    val ksv2 = ("key:subkey", "value2", currentTime)
    val kisv3 = ("key[index]:subkey", "value3", currentTime)
    val kiv4 = ("key[index]", "value4", currentTime)

    val t = Table(
      ("list", "res"),
      (List(kv),  """"key": "value""""),
      (List(kv, ksv2),  """"key": { "subkey": "value2" }"""),
      (List(kv, ksv2, kisv3),  """"key": [ { "subkey": "value3" } ]"""),
      (List(kv, ksv2, kisv3, kiv4),  """"key": [ "value4" ]""")
    )

    forAll(t) { (l, r) => assertMetadataKeyStructure(l, r) }
  }

  it should "return workflow status" in {
    val workflow = WorkflowId.randomId()
    val query = MetadataQuery(workflow, None, Option("status"))
    val events = List(
      MetadataEvent(MetadataKey(workflow, None, "status"), MetadataValue("Running"), Timestamp.valueOf("2000-01-02 12:00:00")),
      MetadataEvent(MetadataKey(workflow, None, "status"), MetadataValue("Done"), Timestamp.valueOf("2000-01-03 12:00:00")),
      MetadataEvent(MetadataKey(workflow, None, "status"), MetadataValue("Submitted"), Timestamp.valueOf("2000-01-01 12:00:00"))
    )
    val expectedRes =
      s"""{
         |"id": "$workflow",
         | "status": "Done"
         |}""".stripMargin

    assertMetadataResponse(GetMetadataQueryAction(query), query, events, expectedRes)
  }

  override def afterAll() = {
    system.shutdown()
  }
}
