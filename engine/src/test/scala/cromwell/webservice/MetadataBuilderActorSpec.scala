package cromwell.webservice

import akka.actor.Props

import java.time.{OffsetDateTime, ZoneOffset}
import akka.pattern.ask
import akka.testkit._
import akka.util.Timeout
import cats.implicits.catsSyntaxValidatedId
import cromwell.core._
import cromwell.services._
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import cromwell.services.metadata.impl.ReadDatabaseMetadataWorkerActor
import cromwell.services.metadata.impl.builder.MetadataBuilderActor
import cromwell.util.AkkaTestUtil.EnhancedTestProbe
import cromwell.webservice.MetadataBuilderActorSpec._
import org.scalatest.flatspec.AsyncFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{Assertion, Succeeded}
import spray.json._

import java.util.UUID
import scala.BigDecimal
import scala.concurrent.Future
import scala.concurrent.duration._
import scala.util.Random

class MetadataBuilderActorSpec
    extends TestKitSuite
    with AsyncFlatSpecLike
    with Matchers
    with TableDrivenPropertyChecks
    with ImplicitSender {

  behavior of "MetadataBuilderActor"

  val defaultSafetyRowNumberThreshold = 1000000
  val defaultTimeout: FiniteDuration = 5.second.dilated
  implicit val timeout: Timeout = defaultTimeout

  def assertMetadataResponse(action: MetadataServiceAction,
                             queryReply: MetadataQuery,
                             events: Seq[MetadataEvent],
                             expectedRes: String,
                             metadataBuilderActorName: String,
                             failedTasks: Boolean = false
  ): Future[Assertion] = {
    val mockReadMetadataWorkerActor = TestProbe("mockReadMetadataWorkerActor")
    def readMetadataWorkerMaker = () => mockReadMetadataWorkerActor.props

    val mba = system.actorOf(
      props = MetadataBuilderActor.props(readMetadataWorkerMaker, 1000000),
      name = metadataBuilderActorName
    )

    val response = mba.ask(action).mapTo[MetadataJsonResponse]
    mockReadMetadataWorkerActor.expectMsg(defaultTimeout, action)
    mockReadMetadataWorkerActor.reply(
      if (failedTasks) FetchFailedJobsMetadataLookupResponse(events) else MetadataLookupResponse(queryReply, events)
    )
    response map { r => r shouldBe a[SuccessfulMetadataJsonResponse] }
    response.mapTo[SuccessfulMetadataJsonResponse] map { b => b.responseJson shouldBe expectedRes.parseJson }
  }

  def assertMetadataFailureResponse(action: MetadataServiceAction,
                                    metadataServiceResponse: MetadataServiceResponse,
                                    expectedException: Exception,
                                    metadataBuilderActorName: String
  ): Future[Assertion] = {
    val mockReadMetadataWorkerActor = TestProbe("mockReadMetadataWorkerActor")
    val mba = system.actorOf(
      props = MetadataBuilderActor.props(() => mockReadMetadataWorkerActor.props, defaultSafetyRowNumberThreshold),
      name = metadataBuilderActorName
    )
    val response = mba.ask(action).mapTo[MetadataServiceResponse]

    mockReadMetadataWorkerActor.expectMsg(defaultTimeout, action)
    mockReadMetadataWorkerActor.reply(metadataServiceResponse)

    response map { r => r shouldBe a[FailedMetadataJsonResponse] }
    response.mapTo[FailedMetadataJsonResponse] map { b =>
      b.reason.getClass shouldBe expectedException.getClass
      b.reason.getMessage shouldBe expectedException.getMessage
    }
  }

  it should "build workflow scope tree from metadata events" in {
    def makeEvent(workflow: WorkflowId, key: Option[MetadataJobKey]) =
      MetadataEvent(MetadataKey(workflow, key, "NOT_CHECKED"), MetadataValue("NOT_CHECKED"))

    val workflowA = WorkflowId.randomId()

    val workflowACalls = List(
      Option(MetadataJobKey("callB", Option(1), 3)),
      Option(MetadataJobKey("callB", None, 1)),
      Option(MetadataJobKey("callB", Option(1), 2)),
      Option(MetadataJobKey("callA", None, 1)),
      Option(MetadataJobKey("callB", Option(1), 1)),
      Option(MetadataJobKey("callB", Option(0), 1)),
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
    assertMetadataResponse(
      action = queryAction,
      queryReply = mdQuery,
      events = workflowAEvents,
      expectedRes = expectedRes,
      metadataBuilderActorName = "mba-scope-tree"
    )
  }

  type EventBuilder = (String, String, OffsetDateTime)

  def makeEvent(
    workflow: WorkflowId
  )(key: String, value: MetadataValue, offsetDateTime: OffsetDateTime): MetadataEvent =
    MetadataEvent(MetadataKey(workflow, None, key), Option(value), offsetDateTime)

  def makeCallEvent(
    workflow: WorkflowId
  )(key: String, value: MetadataValue, offsetDateTime: OffsetDateTime): MetadataEvent = {
    val jobKey = MetadataJobKey("fqn", None, 1)
    MetadataEvent(MetadataKey(workflow, Option(jobKey), key), Option(value), offsetDateTime)
  }

  // noinspection ScalaUnusedSymbol
  def makeEmptyValue(
    workflow: WorkflowId
  )(key: String, value: MetadataValue, offsetDateTime: OffsetDateTime): MetadataEvent =
    MetadataEvent(MetadataKey(workflow, None, key), None, offsetDateTime)

  def assertMetadataKeyStructure(eventList: List[EventBuilder],
                                 expectedJson: String,
                                 workflow: WorkflowId = WorkflowId.randomId(),
                                 eventMaker: WorkflowId => (String, MetadataValue, OffsetDateTime) => MetadataEvent =
                                   makeEvent,
                                 metadataBuilderActorName: String,
                                 isFailedTaskFetch: Boolean = false
  ): Future[Assertion] = {

    val events = eventList map { e => (e._1, MetadataValue(e._2), e._3) } map Function.tupled(eventMaker(workflow))
    val expectedRes = s"""{ "calls": {}, $expectedJson, "id":"$workflow" }"""

    val mdQuery = MetadataQuery(workflow, None, None, None, None, expandSubWorkflows = false)
    val queryAction = GetSingleWorkflowMetadataAction(workflow, None, None, expandSubWorkflows = false)
    assertMetadataResponse(queryAction, mdQuery, events, expectedRes, metadataBuilderActorName)
  }

  it should "build a basic workflow cost response" in {
    val workflowId = WorkflowId.randomId()
    val workflowState = WorkflowSucceeded
    val events = Seq(
      MetadataEvent(MetadataKey(workflowId, Option(MetadataJobKey("callA", None, 1)), CallMetadataKeys.VmStartTime),
                    MetadataValue("2023-10-30T09:00:00Z")
      ),
      MetadataEvent(MetadataKey(workflowId, Option(MetadataJobKey("callA", None, 1)), CallMetadataKeys.VmEndTime),
                    MetadataValue("2023-10-30T11:00:00Z")
      ),
      MetadataEvent(MetadataKey(workflowId, Option(MetadataJobKey("callA", None, 1)), CallMetadataKeys.VmCostPerHour),
                    MetadataValue(BigDecimal(3.5))
      ),
      MetadataEvent(MetadataKey(workflowId, Option(MetadataJobKey("callB", None, 1)), CallMetadataKeys.VmStartTime),
                    MetadataValue("2023-10-30T11:01:00Z")
      ),
      MetadataEvent(MetadataKey(workflowId, Option(MetadataJobKey("callB", None, 1)), CallMetadataKeys.VmEndTime),
                    MetadataValue("2023-10-30T12:31:00Z")
      ),
      MetadataEvent(MetadataKey(workflowId, Option(MetadataJobKey("callB", None, 1)), CallMetadataKeys.VmCostPerHour),
                    MetadataValue(BigDecimal(0.5))
      )
    )
    val query = MetadataQuery(workflowId, None, None, None, None, expandSubWorkflows = false)

    val expectedRes =
      s"""{
         |"cost": 7.75,
         |"currency": "USD",
         |"id": "${workflowId}",
         |"status": "${workflowState.toString}",
         |"errors": []
         |}""".stripMargin

    val action = GetCost(workflowId)

    val mockReadMetadataWorkerActor = TestProbe("mockReadMetadataWorkerActor")
    def readMetadataWorkerMaker = () => mockReadMetadataWorkerActor.props

    val mba = system.actorOf(
      props = MetadataBuilderActor.props(readMetadataWorkerMaker, 1000000),
      name = "mba-cost-builder"
    )

    val response = mba.ask(action).mapTo[MetadataJsonResponse]
    mockReadMetadataWorkerActor.expectMsg(defaultTimeout, action)
    mockReadMetadataWorkerActor.reply(
      CostResponse(workflowId, workflowState, MetadataLookupResponse(query, events))
    )
    response map { r => r shouldBe a[SuccessfulMetadataJsonResponse] }
    response.mapTo[SuccessfulMetadataJsonResponse] map { b => b.responseJson shouldBe expectedRes.parseJson }
  }

  it should "build a workflow cost response with an error" in {
    val workflowId = WorkflowId.randomId()
    val workflowState = WorkflowSucceeded
    val events = Seq(
      MetadataEvent(MetadataKey(workflowId, Option(MetadataJobKey("callA", None, 1)), CallMetadataKeys.VmStartTime),
                    MetadataValue("2023-10-30T09:00:00Z")
      ),
      MetadataEvent(MetadataKey(workflowId, Option(MetadataJobKey("callA", None, 1)), CallMetadataKeys.VmEndTime),
                    MetadataValue("2023-10-30T11:00:00Z")
      ),
      MetadataEvent(
        MetadataKey(workflowId, Option(MetadataJobKey("callA", None, 1)), CallMetadataKeys.VmCostPerHour),
        MetadataValue(BigDecimal(-1)) // indicates an error when computing vmCostPerHour
      ),
      MetadataEvent(MetadataKey(workflowId, Option(MetadataJobKey("callB", None, 1)), CallMetadataKeys.VmStartTime),
                    MetadataValue("2023-10-30T11:01:00Z")
      ),
      MetadataEvent(MetadataKey(workflowId, Option(MetadataJobKey("callB", None, 1)), CallMetadataKeys.VmEndTime),
                    MetadataValue("2023-10-30T12:31:00Z")
      ),
      MetadataEvent(MetadataKey(workflowId, Option(MetadataJobKey("callB", None, 1)), CallMetadataKeys.VmCostPerHour),
                    MetadataValue(BigDecimal(0.5))
      )
    )
    val query = MetadataQuery(workflowId, None, None, None, None, expandSubWorkflows = false)

    val expectedRes =
      s"""{
         |"cost": 0.75,
         |"currency": "USD",
         |"id": "${workflowId}",
         |"status": "${workflowState.toString}",
         |"errors": ["Couldn't find valid vmCostPerHour for callA.-1.1"]
         |}""".stripMargin

    val action = GetCost(workflowId)

    val mockReadMetadataWorkerActor = TestProbe("mockReadMetadataWorkerActor")
    def readMetadataWorkerMaker = () => mockReadMetadataWorkerActor.props

    val mba = system.actorOf(
      props = MetadataBuilderActor.props(readMetadataWorkerMaker, 1000000),
      name = "mba-cost-builder"
    )

    val response = mba.ask(action).mapTo[MetadataJsonResponse]
    mockReadMetadataWorkerActor.expectMsg(defaultTimeout, action)
    mockReadMetadataWorkerActor.reply(
      CostResponse(workflowId, workflowState, MetadataLookupResponse(query, events))
    )
    response map { r => r shouldBe a[SuccessfulMetadataJsonResponse] }
    response.mapTo[SuccessfulMetadataJsonResponse] map { b => b.responseJson shouldBe expectedRes.parseJson }
  }

  it should "build an error workflow cost response" in {
    val workflowId = WorkflowId.randomId()

    val action = GetCost(workflowId)

    val expectedException = new Exception("Oh nooooo :(")

    val mockReadMetadataWorkerActor = TestProbe("mockReadMetadataWorkerActor")
    def readMetadataWorkerMaker = () => mockReadMetadataWorkerActor.props

    val mba = system.actorOf(
      props = MetadataBuilderActor.props(readMetadataWorkerMaker, 1000000),
      name = "mba-cost-builder"
    )

    val response = mba.ask(action).mapTo[MetadataJsonResponse]
    mockReadMetadataWorkerActor.expectMsg(defaultTimeout, action)
    mockReadMetadataWorkerActor.reply(
      CostFailure(workflowId, expectedException)
    )
    response map { r => r shouldBe a[FailedMetadataJsonResponse] }
    response.mapTo[FailedMetadataJsonResponse] map { b =>
      b.reason.getClass shouldBe expectedException.getClass
      b.reason.getMessage shouldBe expectedException.getMessage
    }
  }

  def costMetadataEvents(wfId: WorkflowId,
                         callName: String,
                         shardIndex: Option[Int] = None,
                         attempt: Int = 1,
                         costPerHour: BigDecimal = BigDecimal(1)
  ) = List(
    MetadataEvent(
      MetadataKey(wfId, Option(MetadataJobKey(callName, shardIndex, attempt)), CallMetadataKeys.VmStartTime),
      MetadataValue("2023-10-30T09:00:00Z")
    ),
    MetadataEvent(MetadataKey(wfId, Option(MetadataJobKey(callName, shardIndex, attempt)), CallMetadataKeys.VmEndTime),
                  MetadataValue("2023-10-30T10:00:00Z")
    ),
    MetadataEvent(
      MetadataKey(wfId, Option(MetadataJobKey(callName, shardIndex, attempt)), CallMetadataKeys.VmCostPerHour),
      MetadataValue(costPerHour)
    )
  )

  it should "build a cost response with subworkflow data" in {
    // hard-coding UUIDs to make test failures easier to diagnose
    val mainWorkflowId = WorkflowId(UUID.fromString("00000000-f76d-4af3-b371-5ba580916729"))
    val subWorkflow1Id = WorkflowId(UUID.fromString("11111111-f76d-4af3-b371-5ba580916729"))
    val subWorkflow2Id = WorkflowId(UUID.fromString("22222222-f76d-4af3-b371-5ba580916729"))
    val subWorkflow1aId = WorkflowId(UUID.fromString("1a1a1a1a-f76d-4af3-b371-5ba580916729"))
    val subWorkflow1bId = WorkflowId(UUID.fromString("1b1b1b1b-f76d-4af3-b371-5ba580916729"))

    val workflowSucceededState = WorkflowSucceeded
    val workflowRunningState = WorkflowRunning

    val mainEvents = List(
      MetadataEvent(MetadataKey(mainWorkflowId, Option(MetadataJobKey("wfMain", None, 1)), "subWorkflowId"),
                    MetadataValue(subWorkflow1Id)
      ),
      MetadataEvent(MetadataKey(mainWorkflowId, Option(MetadataJobKey("wfMain", None, 1)), "subWorkflowId"),
                    MetadataValue(subWorkflow2Id)
      )
    ) ++
      costMetadataEvents(mainWorkflowId, "callA", Option(0), 1) ++
      costMetadataEvents(mainWorkflowId, "callA", Option(1), 1)

    val sub1Events = List(
      MetadataEvent(MetadataKey(subWorkflow1Id, Option(MetadataJobKey("wfSub1", None, 1)), "subWorkflowId"),
                    MetadataValue(subWorkflow1aId)
      ),
      MetadataEvent(MetadataKey(subWorkflow1Id, Option(MetadataJobKey("wfSub1", None, 1)), "subWorkflowId"),
                    MetadataValue(subWorkflow1bId)
      )
    )

    val sub2Events =
      costMetadataEvents(subWorkflow2Id, "call2A") ++
        costMetadataEvents(subWorkflow2Id, "call2B") ++
        costMetadataEvents(subWorkflow2Id, "call2C")
    val sub1aEvents = costMetadataEvents(subWorkflow1aId, "call1aA", costPerHour = BigDecimal(-1))
    val sub1bEvents =
      costMetadataEvents(subWorkflow1bId, "call1bA", attempt = 1) ++
        costMetadataEvents(subWorkflow1bId, "call1bA", attempt = 2)

    val mainQuery = MetadataQuery(mainWorkflowId, None, None, None, None, expandSubWorkflows = true)
    val sub1Query = MetadataQuery(subWorkflow1Id, None, None, None, None, expandSubWorkflows = true)
    val sub2Query = MetadataQuery(subWorkflow2Id, None, None, None, None, expandSubWorkflows = true)
    val sub1aQuery = MetadataQuery(subWorkflow1aId, None, None, None, None, expandSubWorkflows = true)
    val sub1bQuery = MetadataQuery(subWorkflow1bId, None, None, None, None, expandSubWorkflows = true)
    val mainAction = GetCost(mainWorkflowId)

    // Mock ReadDatabaseMetadataWorkerActor to return the expected metadata results for each query.
    // Would normally have done this with expect/reply on a TestProbe, but that required the messages to
    // be sent in a deterministic order, which is not the case here.
    class TestReadDatabaseMetadataWorkerActorForCost extends ReadDatabaseMetadataWorkerActor(defaultTimeout, 1000000) {
      override def receive: Receive = {
        case GetCost(wfId) if wfId == mainWorkflowId =>
          sender() ! CostResponse(mainWorkflowId, workflowRunningState, MetadataLookupResponse(mainQuery, mainEvents))
          ()
        case GetCost(wfId) if wfId == subWorkflow1Id =>
          sender() ! CostResponse(subWorkflow1Id, workflowSucceededState, MetadataLookupResponse(sub1Query, sub1Events))
          ()
        case GetCost(wfId) if wfId == subWorkflow2Id =>
          sender() ! CostResponse(subWorkflow2Id, workflowSucceededState, MetadataLookupResponse(sub2Query, sub2Events))
          ()
        case GetCost(wfId) if wfId == subWorkflow1aId =>
          sender() ! CostResponse(subWorkflow1aId,
                                  workflowSucceededState,
                                  MetadataLookupResponse(sub1aQuery, sub1aEvents)
          )
          ()
        case GetCost(wfId) if wfId == subWorkflow1bId =>
          sender() ! CostResponse(subWorkflow1bId,
                                  workflowSucceededState,
                                  MetadataLookupResponse(sub1bQuery, sub1bEvents)
          )
          ()
        case _ => ()
      }
    }

    val expectedRes =
      s"""{
         |"cost": 7,
         |"currency": "USD",
         |"id": "${mainWorkflowId}",
         |"status": "${workflowRunningState.toString}",
         |"errors": ["Couldn't find valid vmCostPerHour for call1aA.-1.1"]
         |}""".stripMargin

    def readMetadataWorkerMaker = () => Props(new TestReadDatabaseMetadataWorkerActorForCost)
    val mba = system.actorOf(
      props = MetadataBuilderActor.props(readMetadataWorkerMaker, 1000000),
      name = "mba-cost-builder"
    )

    val response = mba.ask(mainAction).mapTo[MetadataJsonResponse]
    response map { r => r shouldBe a[SuccessfulMetadataJsonResponse] }
    response.mapTo[SuccessfulMetadataJsonResponse] map { b => b.responseJson shouldBe expectedRes.parseJson }
  }

  it should "compute cost for calls" in {
    val costJsValRows = Table(
      ("jsInput", "expected"),
      ("""{"attempt": 1, "shardIndex": -1, "vmStartTime": "2023-10-30T04:00:00Z", "vmEndTime": "2023-10-30T05:00:00Z", "vmCostPerHour": 0.0567}""",
       BigDecimal(0.0567).validNel
      ),
      ("""{"attempt": 1, "shardIndex": -1, "vmStartTime": "2023-10-30T04:00:00Z", "vmEndTime": "2023-10-30T04:30:00Z", "vmCostPerHour": 0.0567}""",
       BigDecimal(0.02835).validNel
      ),
      ("""{"attempt": 1, "shardIndex": -1, "vmEndTime": "2023-10-30T05:00:00Z", "vmCostPerHour": 0.0567}""",
       BigDecimal(0).validNel
      ),
      ("""{"attempt": 1, "shardIndex": -1, "vmEndTime": "2023-10-30T05:00:00Z", "vmEndTime": "2023-10-30T04:30:00Z"}""",
       BigDecimal(0).validNel
      ),
      (s"""{"attempt": 1, "shardIndex": -1, "vmStartTime": "2023-10-30 05:00:00Z", "vmCostPerHour": 0.0567}""",
       "Text '2023-10-30 05:00:00Z' could not be parsed at index 10".invalidNel
      ),
      (s"""{"attempt": 1, "shardIndex": -1, "vmStartTime": "2023-10-30T05:00:00Z", "vmEndTime": "2023-10-30T04:30:00Z", "vmCostPerHour": -1}""",
       "Couldn't find valid vmCostPerHour for foo.-1.1".invalidNel
      )
    )

    forAll(costJsValRows) { case (jsInput: String, expected) =>
      val jsVal = jsInput.parseJson
      val retVal = MetadataBuilderActor.computeCost("foo", jsVal)
      retVal shouldBe expected
    }
  }

  it should "compute cost for calls that are still running" in {
    val inputJsVal =
      s"""{
         |"attempt": 1,
         |"shardIndex": -1,
         |"vmStartTime": "${OffsetDateTime.now().minusMinutes(5)}",
         |"vmCostPerHour": 12
         |}""".stripMargin.parseJson

    val retVal = MetadataBuilderActor.computeCost("foo", inputJsVal)

    // Just test that the cost is approximately 1 - we don't know
    // exactly how long the test takes to run
    retVal.getOrElse(BigDecimal(0)) should be > BigDecimal(0.9)
    retVal.getOrElse(BigDecimal(100)) should be < BigDecimal(1.1)
  }

  it should "build the call list for failed tasks when prompted" in {

    def makeEvent(workflow: WorkflowId, key: Option[MetadataJobKey]) =
      MetadataEvent(MetadataKey(workflow, key, "NOT_CHECKED"), MetadataValue("NOT_CHECKED"))

    val workflowA = WorkflowId.randomId()

    val workflowACalls = List(
      Option(MetadataJobKey("callB", Option(1), 3)),
      Option(MetadataJobKey("callB", None, 1)),
      Option(MetadataJobKey("callB", Option(1), 2)),
      Option(MetadataJobKey("callA", None, 1)),
      Option(MetadataJobKey("callB", Option(1), 1)),
      Option(MetadataJobKey("callB", Option(0), 1)),
      None
    )
    val workflowAEvents = workflowACalls map {
      makeEvent(workflowA, _)
    }

    val expectedRes =
      s"""{
         |"${workflowA}": {
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
         |  "NOT_CHECKED": "NOT_CHECKED"
         | }
         |}""".stripMargin

    val mdQuery = MetadataQuery(workflowA, None, None, None, None, expandSubWorkflows = false)
    val queryAction = GetMetadataAction(mdQuery)
    assertMetadataResponse(
      action = queryAction,
      queryReply = mdQuery,
      events = workflowAEvents,
      expectedRes = expectedRes,
      failedTasks = true,
      metadataBuilderActorName = "mba-failed-tasks"
    )
  }

  it should "assume the event list is ordered and keep last event if 2 events have same key" in {
    val eventBuilderList = List(
      ("a", "aLater", OffsetDateTime.parse("2000-01-02T12:00:00Z")),
      ("a", "a", OffsetDateTime.parse("2000-01-01T12:00:00Z"))
    )
    val expectedRes =
      """"a": "a"""".stripMargin

    assertMetadataKeyStructure(
      eventList = eventBuilderList,
      expectedJson = expectedRes,
      metadataBuilderActorName = "mba-same-key"
    )
  }

  it should "use CRDT ordering instead of timestamp for workflow state" in {
    val eventBuilderList = List(
      ("status", "Succeeded", OffsetDateTime.now),
      ("status", "Running", OffsetDateTime.now.plusSeconds(1))
    )
    val expectedRes =
      """"status": "Succeeded"""".stripMargin

    assertMetadataKeyStructure(
      eventList = eventBuilderList,
      expectedJson = expectedRes,
      metadataBuilderActorName = "mba-not-workflow-state"
    )
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

    assertMetadataKeyStructure(
      eventList = eventBuilderList,
      expectedJson = expectedRes,
      workflow = workflowId,
      eventMaker = makeCallEvent,
      metadataBuilderActorName = "mba-not-execution-status"
    )
  }

  it should "use reverse date ordering (oldest first) for event start and stop values" in {
    val eventBuilderList = List(
      ("start", "1990-12-20T12:30:00.000Z", OffsetDateTime.now),
      ("start", "1990-12-20T12:30:01.000Z", OffsetDateTime.now.plusSeconds(1)),
      ("end", "2018-06-02T12:30:00.000Z", OffsetDateTime.now.plusSeconds(2)),
      ("end", "2018-06-02T12:30:01.000Z", OffsetDateTime.now.plusSeconds(3))
    )
    val workflowId = WorkflowId.randomId()
    val expectedRes =
      s""""calls": {
         |    "fqn": [{
         |      "attempt": 1,
         |      "end": "2018-06-02T12:30:00.000Z",
         |      "start": "1990-12-20T12:30:00.000Z",
         |      "shardIndex": -1
         |    }]
         |  },
         |  "id": "$workflowId"""".stripMargin

    assertMetadataKeyStructure(
      eventList = eventBuilderList,
      expectedJson = expectedRes,
      workflow = workflowId,
      eventMaker = makeCallEvent,
      metadataBuilderActorName = "mba-start-end-values"
    )
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

    assertMetadataKeyStructure(
      eventList = eventBuilderList,
      expectedJson = expectedRes,
      metadataBuilderActorName = "mba-object-key"
    )
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

    assertMetadataKeyStructure(
      eventList = eventBuilderList,
      expectedJson = expectedRes,
      metadataBuilderActorName = "mba-list-key"
    )
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

    assertMetadataKeyStructure(
      eventList = eventBuilderList,
      expectedJson = expectedRes,
      metadataBuilderActorName = "mba-same-index"
    )
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

    assertMetadataKeyStructure(
      eventList = eventBuilderList,
      expectedJson = expectedRes,
      metadataBuilderActorName = "mba-nest-objects"
    )
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

    assertMetadataKeyStructure(
      eventList = eventBuilderList,
      expectedJson = expectedRes,
      metadataBuilderActorName = "mba-nest-lists"
    )
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

    assertMetadataKeyStructure(
      eventList = eventBuilderList,
      expectedJson = expectedRes,
      eventMaker = makeEmptyValue,
      metadataBuilderActorName = "mba-nest-empty"
    )
  }

  it should "override json values if they can't be merged" in {
    val kv = ("key", "value", OffsetDateTime.now)
    val ksv2 = ("key:subkey", "value2", OffsetDateTime.now.plusSeconds(1))
    val kisv3 = ("key[0]:subkey", "value3", OffsetDateTime.now.plusSeconds(2))
    val kiv4 = ("key[0]", "value4", OffsetDateTime.now.plusSeconds(3))

    val tuples = List(
      ("mba-json-1", List(kv), """"key": "value""""),
      ("mba-json-2", List(kv, ksv2), """"key": { "subkey": "value2" }"""),
      ("mba-json-3", List(kv, ksv2, kisv3), """"key": [ { "subkey": "value3" } ]"""),
      ("mba-json-4", List(kv, ksv2, kisv3, kiv4), """"key": [ "value4" ]""")
    )

    Future.sequence(tuples map { case (metadataBuilderActorName, eventList, expectedJson) =>
      assertMetadataKeyStructure(
        eventList = eventList,
        expectedJson = expectedJson,
        metadataBuilderActorName = metadataBuilderActorName
      )
    }) map { assertions =>
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
    assertMetadataResponse(
      action = queryAction,
      queryReply = mdQuery,
      events = events,
      expectedRes = expectedResponse,
      metadataBuilderActorName = "mba-coerce-type"
    )
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
    assertMetadataResponse(
      action = queryAction,
      queryReply = mdQuery,
      events = events,
      expectedRes = expectedResponse,
      metadataBuilderActorName = "mba-unknown-type"
    )
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
    assertMetadataResponse(
      action = queryAction,
      queryReply = mdQuery,
      events = events,
      expectedRes = expectedResponse,
      metadataBuilderActorName = "mba-coerce-fails"
    )
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
    assertMetadataResponse(
      action = queryAction,
      queryReply = mdQuery,
      events = emptyEvents,
      expectedRes = expectedEmptyResponse,
      metadataBuilderActorName = "mba-empty-values"
    )

    val expectedNonEmptyResponse =
      s"""{
         | "calls": {},
         | "hey": "something",
         | "emptyList": ["something", "something"],
         | "id":"$workflowId"
         |}
      """.stripMargin

    assertMetadataResponse(
      action = queryAction,
      queryReply = mdQuery,
      events = valueEvents,
      expectedRes = expectedNonEmptyResponse,
      metadataBuilderActorName = "mba-non-empty-values"
    )
  }

  it should "expand sub workflow metadata when asked for" in {
    val mainWorkflowId = WorkflowId.randomId()
    val subWorkflowId = WorkflowId.randomId()

    val mainEvents = List(
      MetadataEvent(MetadataKey(mainWorkflowId, Option(MetadataJobKey("callA", None, 1)), "subWorkflowId"),
                    MetadataValue(subWorkflowId)
      )
    )

    val subEvents = List(
      MetadataEvent(MetadataKey(mainWorkflowId, None, "some"), MetadataValue("sub workflow info"))
    )

    val mainQuery = MetadataQuery(mainWorkflowId, None, None, None, None, expandSubWorkflows = true)
    val mainQueryAction = GetMetadataAction(mainQuery)

    val subQuery = MetadataQuery(subWorkflowId, None, None, None, None, expandSubWorkflows = true)
    val subQueryAction = GetMetadataAction(subQuery, checkTotalMetadataRowNumberBeforeQuerying = false)

    val parentProbe = TestProbe("parentProbe")

    val mockReadMetadataWorkerActor = TestProbe("mockReadMetadataWorkerActor")
    def readMetadataWorkerMaker = () => mockReadMetadataWorkerActor.props

    val metadataBuilder =
      TestActorRef(
        props = MetadataBuilderActor.props(readMetadataWorkerMaker, 1000000),
        supervisor = parentProbe.ref,
        name = s"MetadataActor-$mainWorkflowId"
      )
    val response = metadataBuilder.ask(mainQueryAction).mapTo[MetadataJsonResponse]
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

    response map { r => r shouldBe a[SuccessfulMetadataJsonResponse] }
    val bmr = response.mapTo[SuccessfulMetadataJsonResponse]
    bmr map { b => b.responseJson shouldBe expandedRes.parseJson }
  }

  it should "NOT expand sub workflow metadata when NOT asked for" in {
    val mainWorkflowId = WorkflowId.randomId()
    val subWorkflowId = WorkflowId.randomId()

    val mainEvents = List(
      MetadataEvent(MetadataKey(mainWorkflowId, Option(MetadataJobKey("callA", None, 1)), "subWorkflowId"),
                    MetadataValue(subWorkflowId)
      )
    )

    val queryNoExpand = MetadataQuery(mainWorkflowId, None, None, None, None, expandSubWorkflows = false)
    val queryNoExpandAction = GetMetadataAction(queryNoExpand)

    val parentProbe = TestProbe("parentProbe")

    val mockReadMetadataWorkerActor = TestProbe("mockReadMetadataWorkerActor")
    def readMetadataWorkerMaker = () => mockReadMetadataWorkerActor.props

    val metadataBuilder = TestActorRef(
      props = MetadataBuilderActor.props(readMetadataWorkerMaker, 1000000),
      supervisor = parentProbe.ref,
      name = s"MetadataActor-$mainWorkflowId"
    )
    val response = metadataBuilder.ask(queryNoExpandAction).mapTo[MetadataJsonResponse]
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

    response map { r => r shouldBe a[SuccessfulMetadataJsonResponse] }
    val bmr = response.mapTo[SuccessfulMetadataJsonResponse]
    bmr map { b => b.responseJson shouldBe nonExpandedRes.parseJson }

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

    def filterEventsByCall(events: Iterable[MetadataEvent])(call: Call): Iterable[MetadataEvent] =
      events collect {
        case e @ MetadataEvent(MetadataKey(_, Some(MetadataJobKey(n, _, _)), _), _, _) if call.name == n => e
      }

    val calls = List(Foo, Bar, Baz, Qux, Quux, Corge)
    val actuals = calls map filterEventsByCall(actual)
    val expecteds = calls map filterEventsByCall(expectations)

    val matchesExpectations = (actuals zip expecteds) map { case (as, es) =>
      as.toList.map(_.toString).sorted == es.toList.map(_.toString).sorted
    }
    matchesExpectations.reduceLeft(_ && _) shouldBe true
  }

  it should "correctly order statuses (even unused ones)" in {
    val workflowId = WorkflowId.randomId()

    def statusEvent(callName: String, status: String) =
      MetadataEvent(MetadataKey(workflowId, Option(MetadataJobKey(callName, None, 1)), "executionStatus"),
                    MetadataValue(status)
      )

    // Combines standard "setup" statuses plus the conclusion status(es).
    def setupStatusesPlusConclusion(callName: String, conclusionStatuses: String*): Vector[MetadataEvent] = Vector(
      statusEvent(callName, "NotStarted"),
      statusEvent(callName, "WaitingForQueueSpace"),
      statusEvent(callName, "QueuedInCromwell"),
      statusEvent(callName, "Starting"),
      statusEvent(callName, "Running")
    ) ++ conclusionStatuses.map(statusEvent(callName, _))

    /** WARNING!
      * Think twice before removing any of these entries! Even if a status is no longer used, the database can
      * (and probably will!) still contain entries specifying that status.
      */

    val events =
      setupStatusesPlusConclusion("Foo", "Done") ++
        setupStatusesPlusConclusion("Bar", "Aborting", "Aborted") ++
        setupStatusesPlusConclusion("Baz", "Failed") ++
        setupStatusesPlusConclusion("Qux", "RetryableFailure") ++
        setupStatusesPlusConclusion("Quux", "Bypassed") ++
        setupStatusesPlusConclusion("Quuux", "Unstartable")

    val expectedRes =
      s"""{
         |  "calls": {
         |    "Foo": [{ "attempt": 1, "executionStatus": "Done", "shardIndex": -1 }],
         |    "Bar": [{ "attempt": 1, "executionStatus": "Aborted", "shardIndex": -1 }],
         |    "Baz": [{ "attempt": 1, "executionStatus": "Failed", "shardIndex": -1 }],
         |    "Qux": [{ "attempt": 1, "executionStatus": "RetryableFailure", "shardIndex": -1 }],
         |    "Quux": [{ "attempt": 1, "executionStatus": "Bypassed", "shardIndex": -1 }],
         |    "Quuux": [{ "attempt": 1, "executionStatus": "Unstartable", "shardIndex": -1 }]
         |  },
         |  "id": "$workflowId"
         |}""".stripMargin

    val mdQuery = MetadataQuery(workflowId, None, None, None, None, expandSubWorkflows = false)
    val queryAction = GetMetadataAction(mdQuery)

    // The result should always be the same regardless of what order the list arrives in (forward, reverse, random):
    assertMetadataResponse(queryAction, mdQuery, events, expectedRes, "mba-statuses-forward")
    assertMetadataResponse(queryAction, mdQuery, events.reverse, expectedRes, "mba-statuses-reverse")
    assertMetadataResponse(queryAction, mdQuery, Random.shuffle(events), expectedRes, "mba-statuses-random")
  }

  it should "politely refuse building metadata JSON if metadata number of rows is too large" in {
    val workflowId = WorkflowId.randomId()

    val mdQuery = MetadataQuery(workflowId, None, None, None, None, expandSubWorkflows = false)
    val action = GetMetadataAction(mdQuery)

    val metadataRowNumber = 100500
    val expectedException =
      new MetadataTooLargeNumberOfRowsException(workflowId, metadataRowNumber, defaultSafetyRowNumberThreshold)
    assertMetadataFailureResponse(
      action = action,
      metadataServiceResponse = MetadataLookupFailedTooLargeResponse(mdQuery, metadataRowNumber),
      expectedException = expectedException,
      metadataBuilderActorName = "mba-too-large"
    )
  }

  it should "politely refuse building metadata JSON if timeout occurs on attempt to read metadata from database" in {
    val workflowId = WorkflowId.randomId()

    val mdQuery = MetadataQuery(workflowId, None, None, None, None, expandSubWorkflows = false)
    val action = GetMetadataAction(mdQuery)

    val expectedException = new MetadataTooLargeTimeoutException(workflowId)
    assertMetadataFailureResponse(
      action = action,
      metadataServiceResponse = MetadataLookupFailedTimeoutResponse(mdQuery),
      expectedException = expectedException,
      metadataBuilderActorName = "mba-read-timeout"
    )
  }
}

object MetadataBuilderActorSpec {

  val y2k: OffsetDateTime = OffsetDateTime.of(2000, 1, 1, 0, 0, 0, 0, ZoneOffset.UTC)

  case class Interval(start: OffsetDateTime, end: OffsetDateTime) {
    def after: Interval = Interval(start = this.end.plusHours(1), end = this.end.plusHours(2))
  }

  val Interval1: Interval = Interval(y2k, y2k.plusHours(1))
  val Interval2: Interval = Interval1.after
  val Interval3: Interval = Interval2.after
  val Interval4: Interval = Interval3.after
  val Interval5: Interval = Interval4.after
  val Interval6: Interval = Interval5.after
  val Interval7: Interval = Interval6.after
  val Interval8: Interval = Interval7.after
  val Interval9: Interval = Interval8.after

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

  def executionEventKey(workflowId: WorkflowId, call: Call, eventIndex: Int, attr: Attr): MetadataKey =
    MetadataKey(workflowId, Option(MetadataJobKey(call.name, None, 1)), executionEventName(eventIndex, attr))

  def ee(workflowId: WorkflowId, call: Call, eventIndex: Int, attr: Attr, value: Any): MetadataEvent = {
    val metadataValue = value match {
      case g: Grouping => g.name
      case o => o
    }
    new MetadataEvent(executionEventKey(workflowId, call, eventIndex, attr), Option(MetadataValue(metadataValue)), y2k)
  }

  def eventKey(workflowId: WorkflowId, call: Call, attr: Attr): MetadataKey =
    MetadataKey(workflowId, Option(MetadataJobKey(call.name, None, 1)), attr.name)

  def pe(workflowId: WorkflowId, call: Call, attr: Attr, value: Any): MetadataEvent =
    new MetadataEvent(eventKey(workflowId, call, attr), Option(MetadataValue(value)), y2k)
}
