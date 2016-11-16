package cromwell.webservice

import java.time.OffsetDateTime
import java.util.UUID

import akka.testkit._
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import cromwell.webservice.PerRequest.RequestComplete
import cromwell.webservice.metadata.MetadataBuilderActor
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpecLike, Matchers}
import spray.http.{StatusCode, StatusCodes}
import spray.json._

import scala.concurrent.duration._
import scala.language.postfixOps

class MetadataBuilderActorSpec extends TestKitSuite("Metadata") with FlatSpecLike with Matchers
  with TableDrivenPropertyChecks with ImplicitSender {

  behavior of "MetadataParser"

  val defaultTimeout = 200 millis
  val mockServiceRegistry = TestProbe()

  def assertMetadataResponse(action: MetadataServiceAction,
                             queryReply: MetadataQuery,
                             events: Seq[MetadataEvent],
                             expectedRes: String) = {
    val parentProbe = TestProbe()
    val metadataBuilder = TestActorRef(MetadataBuilderActor.props(mockServiceRegistry.ref), parentProbe.ref, s"MetadataActor-${UUID.randomUUID()}")
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

    val workflowACalls = List(
      Option(MetadataJobKey("callB", Some(1), 3)),
      Option(MetadataJobKey("callB", None, 1)),
      Option(MetadataJobKey("callB", Some(1), 2)),
      Option(MetadataJobKey("callA", None, 1)),
      Option(MetadataJobKey("callB", Some(1), 1)),
      Option(MetadataJobKey("callB", Some(0), 1)),
      None
    )
    val workflowAEvents = workflowACalls map { makeEvent(workflowA, _) }

    // We'll use a Query instead of a SingleWorkflowMetadataGet, so we expect the WorkflowID this time:
    val expectedRes =
      s"""{
        |  "calls": {
        |    "callB": [{
        |      "attempt": 1,
        |      "NOT_CHECKED": "NOT_CHECKED",
        |      "shardIndex": -1
        |    }, {
        |      "attempt": 1,
        |      "NOT_CHECKED": "NOT_CHECKED",
        |      "shardIndex": 0
        |    }, {
        |      "attempt": 1,
        |      "NOT_CHECKED": "NOT_CHECKED",
        |      "shardIndex": 1
        |    }, {
        |      "attempt": 2,
        |      "NOT_CHECKED": "NOT_CHECKED",
        |      "shardIndex": 1
        |    }, {
        |      "attempt": 3,
        |      "NOT_CHECKED": "NOT_CHECKED",
        |      "shardIndex": 1
        |    }],
        |    "callA": [{
        |      "attempt": 1,
        |      "NOT_CHECKED": "NOT_CHECKED",
        |      "shardIndex": -1
        |    }]
        |  },
        |  "NOT_CHECKED": "NOT_CHECKED",
        |  "id": "$workflowA"
      |}""".stripMargin

    val mdQuery = MetadataQuery(workflowA, None, None, None, None, expandSubWorkflows = false)
    val queryAction = GetMetadataQueryAction(mdQuery)
    assertMetadataResponse(queryAction, mdQuery, workflowAEvents, expectedRes)
  }

  type EventBuilder = (String, String, OffsetDateTime)

  def makeEvent(workflow: WorkflowId)(key: String, value: MetadataValue, offsetDateTime: OffsetDateTime) = {
    MetadataEvent(MetadataKey(workflow, None, key), Option(value), offsetDateTime)
  }

  def assertMetadataKeyStructure(eventList: List[EventBuilder], expectedJson: String) = {
    val workflow = WorkflowId.randomId()

    val events = eventList map { e => (e._1, MetadataValue(e._2), e._3) } map Function.tupled(makeEvent(workflow))
    val expectedRes = s"""{ "calls": {}, $expectedJson, "id":"$workflow" }"""

    val mdQuery = MetadataQuery(workflow, None, None, None, None, expandSubWorkflows = false)
    val queryAction = GetSingleWorkflowMetadataAction(workflow, None, None, expandSubWorkflows = false)
    assertMetadataResponse(queryAction, mdQuery, events, expectedRes)
  }

  it should "keep event with later timestamp for the same key in metadata" in {
    val eventBuilderList = List(
      ("a", "aLater", OffsetDateTime.parse("2000-01-02T12:00:00Z")),
      ("a", "a", OffsetDateTime.parse("2000-01-01T12:00:00Z"))
    )
    val expectedRes =
      """"a": "aLater"""".stripMargin

    assertMetadataKeyStructure(eventBuilderList, expectedRes)
  }

  it should "use CRDT ordering instead of timestamp for status" in {
    val eventBuilderList = List(
      ("status", "Succeeded", OffsetDateTime.now),
      ("status", "Running", OffsetDateTime.now.plusSeconds(1))
    )
    val expectedRes =
      """"status": "Succeeded"""".stripMargin

    assertMetadataKeyStructure(eventBuilderList, expectedRes)
  }

  it should "build JSON object structure from dotted key syntax" in {
    val eventBuilderList = List(
      ("a:b:c", "abc", OffsetDateTime.now),
      ("b:a", "ba", OffsetDateTime.now),
      ("c", "c", OffsetDateTime.now)
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


  it should "build numerically sorted JSON list structure from dotted key syntax" in {
    val eventBuilderList = List(
      ("l[1]", "l1", OffsetDateTime.now),
      ("l[8]", "l8", OffsetDateTime.now),
      ("l[3]", "l3", OffsetDateTime.now),
      ("l[4]", "l4", OffsetDateTime.now),
      ("l[10]", "l10", OffsetDateTime.now),
      ("l[49]", "l49", OffsetDateTime.now)
    )

    val expectedRes =
      """"l": [
        |    "l1", "l3", "l4", "l8", "l10", "l49"
        |  ]""".stripMargin

    assertMetadataKeyStructure(eventBuilderList, expectedRes)
  }

  it should "override elements with same index in a list if they can't be merged together" in {
    val eventBuilderList = List(
      ("l1[1]", "a", OffsetDateTime.now),
      ("l1[2]", "a", OffsetDateTime.now.plusSeconds(1)),
      ("l1[3]", "a", OffsetDateTime.now.plusSeconds(2)),
      ("l1[2]", "b", OffsetDateTime.now.plusSeconds(3)),
      ("l1[3]", "b", OffsetDateTime.now.plusSeconds(4)),
      ("l1[3]", "c", OffsetDateTime.now.plusSeconds(5))
    )

    val expectedRes =
      """"l1": [
        |    "a", "b", "c"
        |  ]""".stripMargin

    assertMetadataKeyStructure(eventBuilderList, expectedRes)
  }

  it should "nest lists and objects together and respect ordering" in {
    val eventBuilderList = List(
      ("l1[0]:l11[0]:l111[2]", "l10l110l1112", OffsetDateTime.now),
      ("l1[0]:l11[1]:l112[1]", "l10l110l1121", OffsetDateTime.now),
      ("l1[0]:l11[1]:l112[0]", "l10l110l1120", OffsetDateTime.now),
      ("l1[1]:a:l12[2]:b", "l11al122b", OffsetDateTime.now),
      ("l1[0]:l11[0]:l111[0]", "l10l110l1110", OffsetDateTime.now),
      ("l1[0]:l11[0]:l111[1]", "l10l110l1111", OffsetDateTime.now),
      ("l1[0]:l11[2]:l121[0]:a:b", "l10l111l1210ab", OffsetDateTime.now),
      ("l1[1]:a:l12[0]:b", "l11al120b", OffsetDateTime.now),
      ("l1[0]:l11[1]:l112[2]", "l10l110l1122", OffsetDateTime.now),
      ("l1[1]:a:l12[1]", "l11al121", OffsetDateTime.now),
      ("l1[0]:l11[2]:l121[0]:a:c", "l10l111l1210ac", OffsetDateTime.now)
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

  it should "support nested lists" in {
    val eventBuilderList = List(
      ("l[0][0]", "l00", OffsetDateTime.now),
      ("l[0][1]", "l01", OffsetDateTime.now),
      ("l[1][0]:a", "l10a", OffsetDateTime.now),
      ("l[1][1]:b", "l11b", OffsetDateTime.now)
    )

    val expectedRes =
      """"l": [
        |       [
        |        "l00", "l01"
        |       ],
        |       [
        |        { "a": "l10a" }, { "b": "l11b" }
        |       ]
        |     ]""".stripMargin

    assertMetadataKeyStructure(eventBuilderList, expectedRes)
  }

  it should "override json values if they can't be merged" in {
    val kv = ("key", "value", OffsetDateTime.now)
    val ksv2 = ("key:subkey", "value2", OffsetDateTime.now.plusSeconds(1))
    val kisv3 = ("key[0]:subkey", "value3", OffsetDateTime.now.plusSeconds(2))
    val kiv4 = ("key[0]", "value4", OffsetDateTime.now.plusSeconds(3))

    val t = Table(
      ("list", "res"),
      (List(kv),  """"key": "value""""),
      (List(kv, ksv2),  """"key": { "subkey": "value2" }"""),
      (List(kv, ksv2, kisv3),  """"key": [ { "subkey": "value3" } ]"""),
      (List(kv, ksv2, kisv3, kiv4),  """"key": [ "value4" ]""")
    )

    forAll(t) { (l, r) => assertMetadataKeyStructure(l, r) }
  }

  it should "coerce values to supported types" in {
    val workflowId = WorkflowId.randomId()
    val events = List(
      makeEvent(workflowId)("a", MetadataValue(2), OffsetDateTime.now),
      makeEvent(workflowId)("b", MetadataValue(2), OffsetDateTime.now),
      makeEvent(workflowId)("c", MetadataValue(2), OffsetDateTime.now),
      makeEvent(workflowId)("d", MetadataValue(2.9), OffsetDateTime.now),
      makeEvent(workflowId)("e", MetadataValue(2.9), OffsetDateTime.now),
      makeEvent(workflowId)("f", MetadataValue(true), OffsetDateTime.now),
      makeEvent(workflowId)("g", MetadataValue(false), OffsetDateTime.now),
      makeEvent(workflowId)("h", MetadataValue("false"), OffsetDateTime.now)
    )

    val expectedResponse =
      s"""{
          | "calls": {},
          | "a": 2,
          | "b": 2,
          | "c": 2,
          | "d": 2.9,
          | "e": 2.9,
          | "f": true,
          | "g": false,
          | "h": "false",
          | "id":"$workflowId"
          | }
      """.stripMargin

    val mdQuery = MetadataQuery(workflowId, None, None, None, None, expandSubWorkflows = false)
    val queryAction = GetMetadataQueryAction(mdQuery)
    assertMetadataResponse(queryAction, mdQuery, events, expectedResponse)
  }

  it should "fall back to string if the type is unknown" in {
    val workflowId = WorkflowId.randomId()
    case class UnknownClass(v: Int)

    val events = List(
      makeEvent(workflowId)("i", MetadataValue(UnknownClass(50)), OffsetDateTime.now)
    )

    val expectedResponse =
      s"""{
          | "calls": {},
          | "i": "UnknownClass(50)",
          | "id":"$workflowId"
          |}
      """.stripMargin

    val mdQuery = MetadataQuery(workflowId, None, None, None, None, expandSubWorkflows = false)
    val queryAction = GetMetadataQueryAction(mdQuery)
    assertMetadataResponse(queryAction, mdQuery, events, expectedResponse)
  }

  it should "fall back to string if the coercion fails" in {
    val workflowId = WorkflowId.randomId()
    val value = MetadataValue("notAnInt", MetadataInt)
    val events = List(
      makeEvent(workflowId)("i", value, OffsetDateTime.now)
    )

    val expectedResponse =
      s"""{
          | "calls": {},
          | "i": "notAnInt",
          | "id":"$workflowId"
          |}
      """.stripMargin

    val mdQuery = MetadataQuery(workflowId, None, None, None, None, expandSubWorkflows = false)
    val queryAction = GetMetadataQueryAction(mdQuery)
    assertMetadataResponse(queryAction, mdQuery, events, expectedResponse)
  }

  it should "render empty Json" in {
    val workflowId = WorkflowId.randomId()
    val mdQuery = MetadataQuery(workflowId, None, None, None, None, expandSubWorkflows = false)
    val queryAction = GetMetadataQueryAction(mdQuery)
    val expectedEmptyResponse = """{}"""
    assertMetadataResponse(queryAction, mdQuery, List.empty, expectedEmptyResponse)
  }

  it should "render empty values" in {
    val workflowId = WorkflowId.randomId()
    val value = MetadataValue("something")
    val emptyEvents = List(
      MetadataEvent.empty(MetadataKey(workflowId, None, "hey")),
      MetadataEvent.empty(MetadataKey(workflowId, None, "emptyList[]"))
    )
    val valueEvents = List(
      MetadataEvent.empty(MetadataKey(workflowId, None, "hey")),
      MetadataEvent.empty(MetadataKey(workflowId, None, "emptyList[]")),
      MetadataEvent(MetadataKey(workflowId, None, "hey"), Option(value), OffsetDateTime.now().plusSeconds(1L)),
      MetadataEvent(MetadataKey(workflowId, None, "emptyList[0]"), Option(value), OffsetDateTime.now().plusSeconds(1L)),
      MetadataEvent(MetadataKey(workflowId, None, "emptyList[1]"), Option(value), OffsetDateTime.now().plusSeconds(1L))
    )

    val expectedEmptyResponse =
      s"""{
          | "calls": {},
          | "hey": {},
          | "emptyList": [],
          | "id":"$workflowId"
          |}
      """.stripMargin

    val mdQuery = MetadataQuery(workflowId, None, None, None, None, expandSubWorkflows = false)
    val queryAction = GetMetadataQueryAction(mdQuery)
    assertMetadataResponse(queryAction, mdQuery, emptyEvents, expectedEmptyResponse)

    val expectedNonEmptyResponse =
      s"""{
          | "calls": {},
          | "hey": "something",
          | "emptyList": ["something", "something"],
          | "id":"$workflowId"
          |}
      """.stripMargin

    assertMetadataResponse(queryAction, mdQuery, valueEvents, expectedNonEmptyResponse)
  }
  
  it should "expand sub workflow metadata when asked for" in {
    val mainWorkflowId = WorkflowId.randomId()
    val subWorkflowId = WorkflowId.randomId()
    
    val mainEvents = List(
      MetadataEvent(MetadataKey(mainWorkflowId, Option(MetadataJobKey("callA", None, 1)), "subWorkflowId"), MetadataValue(subWorkflowId))
    )

    val subEvents = List(
      MetadataEvent(MetadataKey(mainWorkflowId, None, "some"), MetadataValue("sub workflow info"))
    )
    
    val mainQuery = MetadataQuery(mainWorkflowId, None, None, None, None, expandSubWorkflows = true)
    val mainQueryAction = GetMetadataQueryAction(mainQuery)
    
    val subQuery = MetadataQuery(subWorkflowId, None, None, None, None, expandSubWorkflows = true)
    val subQueryAction = GetMetadataQueryAction(subQuery)
    
    val parentProbe = TestProbe()
    val metadataBuilder = TestActorRef(MetadataBuilderActor.props(mockServiceRegistry.ref), parentProbe.ref, s"MetadataActor-${UUID.randomUUID()}")
    metadataBuilder ! mainQueryAction
    mockServiceRegistry.expectMsg(defaultTimeout, mainQueryAction)
    mockServiceRegistry.reply(MetadataLookupResponse(mainQuery, mainEvents))
    mockServiceRegistry.expectMsg(defaultTimeout, subQueryAction)
    mockServiceRegistry.reply(MetadataLookupResponse(subQuery, subEvents))
    
    val expandedRes =
      s"""
         |{
         |  "calls": {
         |    "callA": [
         |      {
         |        "subWorkflowMetadata": {
         |          "some": "sub workflow info",
         |          "calls": {},
         |          "id": "$subWorkflowId"
         |         },
         |         "attempt": 1,
         |         "shardIndex": -1
         |      }
         |    ]
         |  },
         |  "id": "$mainWorkflowId"
         |}
       """.stripMargin

    parentProbe.expectMsgPF(defaultTimeout) {
      case response: RequestComplete[(StatusCode, JsObject)] @unchecked =>
        response.response._1 shouldBe StatusCodes.OK
        response.response._2 shouldBe expandedRes.parseJson
    }
  }
  
  it should "NOT expand sub workflow metadata when NOT asked for" in {
    val mainWorkflowId = WorkflowId.randomId()
    val subWorkflowId = WorkflowId.randomId()

    val mainEvents = List(
      MetadataEvent(MetadataKey(mainWorkflowId, Option(MetadataJobKey("callA", None, 1)), "subWorkflowId"), MetadataValue(subWorkflowId))
    )

    val queryNoExpand = MetadataQuery(mainWorkflowId, None, None, None, None, expandSubWorkflows = false)
    val queryNoExpandAction = GetMetadataQueryAction(queryNoExpand)
    
    val parentProbe = TestProbe()
    val metadataBuilder = TestActorRef(MetadataBuilderActor.props(mockServiceRegistry.ref), parentProbe.ref, s"MetadataActor-${UUID.randomUUID()}")
    metadataBuilder ! queryNoExpandAction
    mockServiceRegistry.expectMsg(defaultTimeout, queryNoExpandAction)
    mockServiceRegistry.reply(MetadataLookupResponse(queryNoExpand, mainEvents))


    val nonExpandedRes =
      s"""
         |{
         |  "calls": {
         |    "callA": [
         |      {
         |        "subWorkflowId": "$subWorkflowId",
         |        "attempt": 1,
         |        "shardIndex": -1
         |      }  
         |    ]
         |  },
         |  "id": "$mainWorkflowId"
         |}
       """.stripMargin
    
    parentProbe.expectMsgPF(defaultTimeout) {
      case response: RequestComplete[(StatusCode, JsObject)] @unchecked =>
        response.response._1 shouldBe StatusCodes.OK
        response.response._2 shouldBe nonExpandedRes.parseJson
    }
  }
}
