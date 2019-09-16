package cromwell.webservice

import java.time.{OffsetDateTime, ZoneOffset}
import java.util.UUID

import akka.pattern.ask
import akka.testkit._
import akka.util.Timeout
import cromwell.core._
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import cromwell.services.metadata.impl.builder.MetadataBuilderActor
import cromwell.services.metadata.impl.builder.MetadataBuilderActor.{BuiltMetadataResponse, MetadataBuilderActorResponse}
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{Assertion, AsyncFlatSpecLike, Matchers, Succeeded}
import org.specs2.mock.Mockito
import spray.json._
import cromwell.util.AkkaTestUtil.EnhancedTestProbe

import scala.concurrent.Future
import scala.concurrent.duration._
import cromwell.webservice.MetadataBuilderActorSpec._


class MetadataBuilderActorSpec extends TestKitSuite("Metadata") with AsyncFlatSpecLike with Matchers with Mockito
  with TableDrivenPropertyChecks with ImplicitSender {

  behavior of "MetadataBuilderActor"

  val defaultTimeout: FiniteDuration = 1.second.dilated
  implicit val timeout: Timeout = defaultTimeout

  def assertMetadataResponse(action: MetadataServiceAction,
                             queryReply: MetadataQuery,
                             events: Seq[MetadataEvent],
                             expectedRes: String): Future[Assertion] = {
    val mockReadMetadataWorkerActor = TestProbe()
    def readMetadataWorkerMaker = () => mockReadMetadataWorkerActor.props


    val mba = system.actorOf(MetadataBuilderActor.props(readMetadataWorkerMaker))
    val response = mba.ask(action).mapTo[MetadataBuilderActorResponse]
    mockReadMetadataWorkerActor.expectMsg(defaultTimeout, action)
    mockReadMetadataWorkerActor.reply(MetadataLookupResponse(queryReply, events))
    response map { r => r shouldBe a [BuiltMetadataResponse] }
    response.mapTo[BuiltMetadataResponse] map { b => b.responseJson shouldBe expectedRes.parseJson}
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
    val queryAction = GetMetadataAction(mdQuery)
    assertMetadataResponse(queryAction, mdQuery, workflowAEvents, expectedRes)
  }

  type EventBuilder = (String, String, OffsetDateTime)

  def makeEvent(workflow: WorkflowId)(key: String, value: MetadataValue, offsetDateTime: OffsetDateTime): MetadataEvent = {
    MetadataEvent(MetadataKey(workflow, None, key), Option(value), offsetDateTime)
  }

  def makeCallEvent(workflow: WorkflowId)(key: String, value: MetadataValue, offsetDateTime: OffsetDateTime) = {
    val jobKey = MetadataJobKey("fqn", None, 1)
    MetadataEvent(MetadataKey(workflow, Option(jobKey), key), Option(value), offsetDateTime)
  }

  def makeEmptyValue(workflow: WorkflowId)(key: String, value: MetadataValue, offsetDateTime: OffsetDateTime) = {
    MetadataEvent(MetadataKey(workflow, None, key), None, offsetDateTime)
  }

  def assertMetadataKeyStructure(eventList: List[EventBuilder],
                                 expectedJson: String,
                                 workflow: WorkflowId = WorkflowId.randomId(),
                                 eventMaker: WorkflowId => (String, MetadataValue, OffsetDateTime) => MetadataEvent = makeEvent) = {

    val events = eventList map { e => (e._1, MetadataValue(e._2), e._3) } map Function.tupled(eventMaker(workflow))
    val expectedRes = s"""{ "calls": {}, $expectedJson, "id":"$workflow" }"""

    val mdQuery = MetadataQuery(workflow, None, None, None, None, expandSubWorkflows = false)
    val queryAction = GetSingleWorkflowMetadataAction(workflow, None, None, expandSubWorkflows = false)
    assertMetadataResponse(queryAction, mdQuery, events, expectedRes)
  }

  it should "assume the event list is ordered and keep last event if 2 events have same key" in {
    val eventBuilderList = List(
      ("a", "aLater", OffsetDateTime.parse("2000-01-02T12:00:00Z")),
      ("a", "a", OffsetDateTime.parse("2000-01-01T12:00:00Z"))
    )
    val expectedRes =
      """"a": "a"""".stripMargin

    assertMetadataKeyStructure(eventBuilderList, expectedRes)
  }

  it should "use CRDT ordering instead of timestamp for workflow state" in {
    val eventBuilderList = List(
      ("status", "Succeeded", OffsetDateTime.now),
      ("status", "Running", OffsetDateTime.now.plusSeconds(1))
    )
    val expectedRes =
      """"status": "Succeeded"""".stripMargin

    assertMetadataKeyStructure(eventBuilderList, expectedRes)
  }

  it should "use CRDT ordering instead of timestamp for call execution status" in {
    val eventBuilderList = List(
      ("executionStatus", "Done", OffsetDateTime.now),
      ("executionStatus", "Running", OffsetDateTime.now.plusSeconds(1))
    )
    val workflowId = WorkflowId.randomId()
    val expectedRes =
      s""""calls": {
        |    "fqn": [{
        |      "attempt": 1,
        |      "executionStatus": "Done",
        |      "shardIndex": -1
        |    }]
        |  },
        |  "id": "$workflowId"""".stripMargin

    assertMetadataKeyStructure(eventBuilderList, expectedRes, workflowId, makeCallEvent)
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

  it should "support nested empty lists" in {
    val eventBuilderList = List(
      ("l[0][]", null, OffsetDateTime.now),
      ("l[1][]", null, OffsetDateTime.now)
    )

    val expectedRes =
      """"l": [
        |       [], []
        |     ]""".stripMargin

    assertMetadataKeyStructure(eventBuilderList, expectedRes, eventMaker = makeEmptyValue)
  }

  it should "override json values if they can't be merged" in {
    val kv = ("key", "value", OffsetDateTime.now)
    val ksv2 = ("key:subkey", "value2", OffsetDateTime.now.plusSeconds(1))
    val kisv3 = ("key[0]:subkey", "value3", OffsetDateTime.now.plusSeconds(2))
    val kiv4 = ("key[0]", "value4", OffsetDateTime.now.plusSeconds(3))

    val t = List(
      (List(kv),  """"key": "value""""),
      (List(kv, ksv2),  """"key": { "subkey": "value2" }"""),
      (List(kv, ksv2, kisv3),  """"key": [ { "subkey": "value3" } ]"""),
      (List(kv, ksv2, kisv3, kiv4),  """"key": [ "value4" ]""")
    )

    Future.sequence(t map { case (l, r) => assertMetadataKeyStructure(l, r) }) map { assertions =>
      assertions should contain only Succeeded
    }
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
    val queryAction = GetMetadataAction(mdQuery)
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
    val queryAction = GetMetadataAction(mdQuery)
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
    val queryAction = GetMetadataAction(mdQuery)
    assertMetadataResponse(queryAction, mdQuery, events, expectedResponse)
  }

  it should "render empty Json" in {
    val workflowId = WorkflowId.randomId()
    val mdQuery = MetadataQuery(workflowId, None, None, None, None, expandSubWorkflows = false)
    val queryAction = GetMetadataAction(mdQuery)
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
    val queryAction = GetMetadataAction(mdQuery)
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
    val mainQueryAction = GetMetadataAction(mainQuery)
    
    val subQuery = MetadataQuery(subWorkflowId, None, None, None, None, expandSubWorkflows = true)
    val subQueryAction = GetMetadataAction(subQuery)
    
    val parentProbe = TestProbe()

    val mockReadMetadataWorkerActor = TestProbe()
    def readMetadataWorkerMaker = () => mockReadMetadataWorkerActor.props

    val metadataBuilder = TestActorRef(MetadataBuilderActor.props(readMetadataWorkerMaker), parentProbe.ref, s"MetadataActor-${UUID.randomUUID()}")
    val response = metadataBuilder.ask(mainQueryAction).mapTo[MetadataBuilderActorResponse]
    mockReadMetadataWorkerActor.expectMsg(defaultTimeout, mainQueryAction)
    mockReadMetadataWorkerActor.reply(MetadataLookupResponse(mainQuery, mainEvents))
    mockReadMetadataWorkerActor.expectMsg(defaultTimeout, subQueryAction)
    mockReadMetadataWorkerActor.reply(MetadataLookupResponse(subQuery, subEvents))
    
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

    response map { r => r shouldBe a [BuiltMetadataResponse] }
    val bmr = response.mapTo[BuiltMetadataResponse]
    bmr map { b => b.responseJson shouldBe expandedRes.parseJson}
  }
  
  it should "NOT expand sub workflow metadata when NOT asked for" in {
    val mainWorkflowId = WorkflowId.randomId()
    val subWorkflowId = WorkflowId.randomId()

    val mainEvents = List(
      MetadataEvent(MetadataKey(mainWorkflowId, Option(MetadataJobKey("callA", None, 1)), "subWorkflowId"), MetadataValue(subWorkflowId))
    )

    val queryNoExpand = MetadataQuery(mainWorkflowId, None, None, None, None, expandSubWorkflows = false)
    val queryNoExpandAction = GetMetadataAction(queryNoExpand)
    
    val parentProbe = TestProbe()

    val mockReadMetadataWorkerActor = TestProbe()
    def readMetadataWorkerMaker= () => mockReadMetadataWorkerActor.props

    val metadataBuilder = TestActorRef(MetadataBuilderActor.props(readMetadataWorkerMaker), parentProbe.ref, s"MetadataActor-${UUID.randomUUID()}")
    val response = metadataBuilder.ask(queryNoExpandAction).mapTo[MetadataBuilderActorResponse]
    mockReadMetadataWorkerActor.expectMsg(defaultTimeout, queryNoExpandAction)
    mockReadMetadataWorkerActor.reply(MetadataLookupResponse(queryNoExpand, mainEvents))


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

    response map { r => r shouldBe a [BuiltMetadataResponse] }
    val bmr = response.mapTo[BuiltMetadataResponse]
    bmr map { b => b.responseJson shouldBe nonExpandedRes.parseJson}

  }

  it should "group execution events properly" in {
    val workflowId = WorkflowId.randomId()

    val events = List(
      // Foo should produce a synthetic grouping event with span [Interval1.Start, Interval2.end].
      ee(workflowId, Foo, eventIndex = 1, StartTime, Interval2.start),
      ee(workflowId, Foo, eventIndex = 1, EndTime, Interval2.end),
      ee(workflowId, Foo, eventIndex = 1, Grouping, Localizing),
      ee(workflowId, Foo, eventIndex = 2, StartTime, Interval1.start),
      ee(workflowId, Foo, eventIndex = 2, EndTime, Interval1.end),
      ee(workflowId, Foo, eventIndex = 2, Grouping, Localizing),

      // Bar should produce a synthetic grouping event with span [Interval3.Start, None] (no closing end timestamp for the last event)
      ee(workflowId, Bar, eventIndex = 3, StartTime, Interval5.start),
      ee(workflowId, Bar, eventIndex = 3, Grouping, Delocalizing),
      ee(workflowId, Bar, eventIndex = 4, StartTime, Interval3.start),
      ee(workflowId, Bar, eventIndex = 4, EndTime, Interval3.end),
      ee(workflowId, Bar, eventIndex = 4, Grouping, Delocalizing),
      ee(workflowId, Bar, eventIndex = 5, StartTime, Interval4.start),
      ee(workflowId, Bar, eventIndex = 5, EndTime, Interval4.end),
      ee(workflowId, Bar, eventIndex = 5, Grouping, Delocalizing),

      ee(workflowId, Baz, eventIndex = 6, StartTime, Interval2.start),
      ee(workflowId, Baz, eventIndex = 6, EndTime, Interval2.end),
      ee(workflowId, Baz, eventIndex = 6, Grouping, Localizing),
      ee(workflowId, Baz, eventIndex = 7, StartTime, Interval1.start),
      ee(workflowId, Baz, eventIndex = 7, EndTime, Interval1.end),
      ee(workflowId, Baz, eventIndex = 7, Grouping, Localizing),

      // Qux should not get grouped.
      ee(workflowId, Qux, eventIndex = 8, StartTime, Interval7.start),
      ee(workflowId, Qux, eventIndex = 8, EndTime, Interval7.end),

      // These 'plain events' don't have groups and aren't executionEvents. These are just here to make sure their presence doesn't
      // cause problems.
      pe(workflowId, Quux, StartTime, Interval8.start),
      pe(workflowId, Quux, EndTime, Interval8.end),

      pe(workflowId, Corge, StartTime, Interval9.start),
      pe(workflowId, Corge, EndTime, Interval9.end)
    )

    // The values for `eventIndex` are purely artifacts of how the logic in `groupEvents` chooses one of the keys
    // in the real execution events to use as the key for the synthetic execution events. If that logic changes
    // these expectations should be updated. Whatever logic is used should probably use one of the keys from the
    // real execution events to prevent colliding with any other set of events.
    val expectations = Set(
      ee(workflowId, Foo, eventIndex = 1, StartTime, Interval1.start),
      ee(workflowId, Foo, eventIndex = 1, EndTime, Interval2.end),
      ee(workflowId, Foo, eventIndex = 1, Description, Localizing.name),

      ee(workflowId, Bar, eventIndex = 3, StartTime, Interval3.start),
      ee(workflowId, Bar, eventIndex = 3, Description, Delocalizing.name),

      ee(workflowId, Baz, eventIndex = 6, StartTime, Interval1.start),
      ee(workflowId, Baz, eventIndex = 6, EndTime, Interval2.end),
      ee(workflowId, Baz, eventIndex = 6, Description, Localizing.name),

      ee(workflowId, Qux, eventIndex = 8, StartTime, Interval7.start),
      ee(workflowId, Qux, eventIndex = 8, EndTime, Interval7.end),

      pe(workflowId, Quux, StartTime, Interval8.start),
      pe(workflowId, Quux, EndTime, Interval8.end),

      pe(workflowId, Corge, StartTime, Interval9.start),
      pe(workflowId, Corge, EndTime, Interval9.end)
    )

    val actual = MetadataBuilderActor.groupEvents(events).toSet

    def filterEventsByCall(events: Iterable[MetadataEvent])(call: Call): Iterable[MetadataEvent] = {
      events collect { case e@MetadataEvent(MetadataKey(_, Some(MetadataJobKey(n, _, _)), _), _, _) if call.name == n => e}
    }

    val calls = List(Foo, Bar, Baz, Qux, Quux, Corge)
    val actuals = calls map filterEventsByCall(actual)
    val expecteds = calls map filterEventsByCall(expectations)

    val matchesExpectations = (actuals zip expecteds) map {
      case (as, es) => (as.toList.map { _.toString } sorted) == (es.toList.map { _.toString } sorted)
    }
    matchesExpectations.reduceLeft(_ && _) shouldBe true
  }
}

object MetadataBuilderActorSpec {

  val y2k = OffsetDateTime.of(2000, 1, 1, 0, 0, 0, 0, ZoneOffset.UTC)

  case class Interval(start: OffsetDateTime, end: OffsetDateTime) {
    def after: Interval = Interval(start = this.end.plusHours(1), end = this.end.plusHours(2))
  }

  val Interval1 = Interval(y2k, y2k.plusHours(1))
  val Interval2 = Interval1.after
  val Interval3 = Interval2.after
  val Interval4 = Interval3.after
  val Interval5 = Interval4.after
  val Interval6 = Interval5.after
  val Interval7 = Interval6.after
  val Interval8 = Interval7.after
  val Interval9 = Interval8.after

  sealed trait Attr {
    val name: String
  }

  case object StartTime extends Attr {
    override val name = "startTime"
  }

  case object EndTime extends Attr {
    override val name = "endTime"
  }

  case object Grouping extends Attr {
    override val name = "grouping"
  }

  case object Description extends Attr {
    override val name = "description"
  }

  sealed trait Call {
    def name: String = getClass.getSimpleName.toLowerCase.takeWhile(_.isLetter)
  }
  case object Foo extends Call
  case object Bar extends Call
  case object Baz extends Call
  case object Qux extends Call
  case object Quux extends Call
  case object Corge extends Call

  sealed trait Grouping {
    def name: String = getClass.getSimpleName.takeWhile(_.isLetter)
  }
  case object Localizing extends Grouping
  case object Delocalizing extends Grouping

  def executionEventName(i: Int, a: Attr): String = s"executionEvents[$i]:${a.name}"

  def executionEventKey(workflowId: WorkflowId, call: Call, eventIndex: Int, attr: Attr): MetadataKey = MetadataKey(workflowId, Option(MetadataJobKey(call.name, None, 1)), executionEventName(eventIndex, attr))


  def ee(workflowId: WorkflowId, call: Call, eventIndex: Int, attr: Attr, value: Any): MetadataEvent = {
    val metadataValue = value match {
      case g: Grouping => g.name
      case o => o
    }
    new MetadataEvent(executionEventKey(workflowId, call, eventIndex, attr), Option(MetadataValue(metadataValue)), y2k)
  }

  def eventKey(workflowId: WorkflowId, call: Call, attr: Attr): MetadataKey = MetadataKey(workflowId, Option(MetadataJobKey(call.name, None, 1)), attr.name)

  def pe(workflowId: WorkflowId, call: Call, attr: Attr, value: Any): MetadataEvent = new MetadataEvent(eventKey(workflowId, call, attr), Option(MetadataValue(value)), y2k)
}
