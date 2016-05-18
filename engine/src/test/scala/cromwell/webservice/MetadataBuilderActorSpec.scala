package cromwell.webservice

import akka.actor.ActorSystem
import akka.testkit._
import cromwell.core.WorkflowId
import cromwell.services.MetadataServiceActor._
import cromwell.webservice.PerRequest.RequestComplete
import org.joda.time.DateTime
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{BeforeAndAfterAll, FlatSpecLike, Matchers}
import spray.http.{StatusCode, StatusCodes}
import spray.json._
import scala.concurrent.duration._
import scala.language.postfixOps

class MetadataBuilderActorSpec extends TestKit(ActorSystem("Metadata")) with FlatSpecLike with Matchers with TableDrivenPropertyChecks with ImplicitSender with BeforeAndAfterAll {

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
      MetadataEvent(MetadataKey(workflow, key, "NOT_CHECKED"), MetadataValue("NOT_CHECKED"))
    }

    val workflowA = WorkflowId.randomId()
    val workflowB = WorkflowId.randomId()

    val workflowACalls = List(
      Option(MetadataJobKey("callB", Some(1), 3)),
      Option(MetadataJobKey("callB", None, 1)),
      Option(MetadataJobKey("callB", Some(1), 2)),
      Option(MetadataJobKey("callA", None, 1)),
      Option(MetadataJobKey("callB", Some(1), 1)),
      Option(MetadataJobKey("callB", Some(0), 1)),
      None
    )
    val workflowBCalls = List(Option(MetadataJobKey("callA", None, 1)), None)

    val workflowAEvents = workflowACalls map { makeEvent(workflowA, _) }
    val workflowBEvents = workflowBCalls map { makeEvent(workflowB, _) }

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
        |  },
        |  "$workflowB": {
        |    "calls": {
        |      "callA": [{
        |        "attempt": 1,
        |        "NOT_CHECKED": "NOT_CHECKED",
        |        "shardIndex": -1
        |      }]
        |    },
        |    "NOT_CHECKED": "NOT_CHECKED"
        |  }
        |}""".stripMargin

    val mdQuery = MetadataQuery(None, None, None)
    val queryAction = GetMetadataQueryAction(mdQuery)
    assertMetadataResponse(queryAction, mdQuery, workflowAEvents ++ workflowBEvents, expectedRes)
  }

  type EventBuilder = (String, String, DateTime)

  def assertMetadataKeyStructure(eventList: List[EventBuilder], expectedJson: String) = {
    val workflow = WorkflowId.randomId()

    def makeEvent(workflow: WorkflowId)(key: String, value: String, timestamp: DateTime) = {
      MetadataEvent(MetadataKey(workflow, None, key), MetadataValue(value), timestamp)
    }

    val events = eventList map Function.tupled(makeEvent(workflow))
    val expectedRes = s"""{ "calls": {}, $expectedJson } """

    val mdQuery = MetadataQuery(Option(workflow), None, None)
    val queryAction = GetAllMetadataAction(workflow)
    assertMetadataResponse(queryAction, mdQuery, events, expectedRes)
  }

  it should "keep event with later timestamp for the same key in metadata" in {
    val eventBuilderList = List(
      ("a", "aLater", DateTime.parse("2000-01-02")),
      ("a", "a", DateTime.parse("2000-01-01"))
    )
    val expectedRes =
      """"a": "aLater"""".stripMargin

    assertMetadataKeyStructure(eventBuilderList, expectedRes)
  }

  it should "build JSON object structure from dotted key syntax" in {
    val eventBuilderList = List(
      ("a:b:c", "abc", DateTime.now),
      ("b:a", "ba", DateTime.now),
      ("c", "c", DateTime.now)
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
      ("l1[2]", "l12", DateTime.now),
      ("l1[k]", "l1k", DateTime.now),
      ("l1[3]", "l13", DateTime.now),
      ("l1[i]", "l1i", DateTime.now),
      ("l1[10]", "l110", DateTime.now),
      ("l1[j]", "l1j", DateTime.now)
    )

    val expectedRes =
      """"l1": [
        |    "l110", "l12", "l13", "l1i", "l1j", "l1k"
        |  ]""".stripMargin

    assertMetadataKeyStructure(eventBuilderList, expectedRes)
  }

  it should "override elements with same index in a list if they can't be merged together" in {
    val eventBuilderList = List(
      ("l1[1]", "a", DateTime.now),
      ("l1[2]", "a", DateTime.now),
      ("l1[3]", "a", DateTime.now),
      ("l1[2]", "b", DateTime.now),
      ("l1[3]", "b", DateTime.now),
      ("l1[3]", "c", DateTime.now)
    )

    val expectedRes =
      """"l1": [
        |    "a", "b", "c"
        |  ]""".stripMargin

    assertMetadataKeyStructure(eventBuilderList, expectedRes)
  }

  it should "nest lists and objects together and respect ordering" in {
    val eventBuilderList = List(
      ("l1[0]:l11[0]:l111[2]", "l10l110l1112", DateTime.now),
      ("l1[0]:l11[1]:l112[1]", "l10l110l1121", DateTime.now),
      ("l1[0]:l11[1]:l112[0]", "l10l110l1120", DateTime.now),
      ("l1[1]:a:l12[2]:b", "l11al122b", DateTime.now),
      ("l1[0]:l11[0]:l111[0]", "l10l110l1110", DateTime.now),
      ("l1[0]:l11[0]:l111[1]", "l10l110l1111", DateTime.now),
      ("l1[0]:l11[2]:l121[0]:a:b", "l10l111l1210ab", DateTime.now),
      ("l1[1]:a:l12[0]:b", "l11al120b", DateTime.now),
      ("l1[0]:l11[1]:l112[2]", "l10l110l1122", DateTime.now),
      ("l1[1]:a:l12[1]", "l11al121", DateTime.now),
      ("l1[0]:l11[2]:l121[0]:a:c", "l10l111l1210ac", DateTime.now)
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
    val kv = ("key", "value", DateTime.now)
    val ksv2 = ("key:subkey", "value2", DateTime.now)
    val kisv3 = ("key[index]:subkey", "value3", DateTime.now)
    val kiv4 = ("key[index]", "value4", DateTime.now)

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
    val query = MetadataQuery(Option(workflow), None, Option("status"))
    val events = List(
      MetadataEvent(MetadataKey(workflow, None, "status"), MetadataValue("Running"), DateTime.parse("2000-01-02")),
      MetadataEvent(MetadataKey(workflow, None, "status"), MetadataValue("Done"), DateTime.parse("2000-01-03")),
      MetadataEvent(MetadataKey(workflow, None, "status"), MetadataValue("Submitted"), DateTime.parse("2000-01-01"))
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
