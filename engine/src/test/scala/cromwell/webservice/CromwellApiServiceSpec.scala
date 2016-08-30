package cromwell.webservice

import java.time.OffsetDateTime
import java.util.UUID

import akka.actor.{Actor, Props}
import akka.pattern.ask
import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestkitSpec
import cromwell.core.Tags._
import cromwell.core._
import cromwell.engine.workflow.WorkflowDescriptorBuilder
import cromwell.engine.workflow.WorkflowManagerActor.AbortWorkflowCommand
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor.{BatchSubmitWorkflows, SubmitWorkflow, WorkflowSubmittedToStore, WorkflowsBatchSubmittedToStore}
import cromwell.server.{CromwellServerActor, CromwellSystem}
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import cromwell.services.metadata.impl.MetadataSummaryRefreshActor.MetadataSummarySuccess
import cromwell.util.SampleWdl.HelloWorld
import cromwell.webservice.CromwellApiHandler._
import org.scalatest.concurrent.{PatienceConfiguration, ScalaFutures}
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito
import spray.http.{DateTime => _, _}
import spray.json.DefaultJsonProtocol._
import spray.json._
import spray.testkit.ScalatestRouteTest

object MockWorkflowStoreActor {
  val submittedWorkflowId = WorkflowId(UUID.randomUUID())
}

class MockWorkflowStoreActor extends Actor {
  import MockWorkflowStoreActor.submittedWorkflowId

  override def receive = {
    case SubmitWorkflow(source) => sender ! WorkflowSubmittedToStore(submittedWorkflowId)
    case BatchSubmitWorkflows(sources) =>
      val response = WorkflowsBatchSubmittedToStore(sources map { _ => submittedWorkflowId })
      sender ! response
  }
}

object MockWorkflowManagerActor {
  val runningWorkflowId = WorkflowId.randomId()
  val unknownId = WorkflowId.randomId()
  val submittedScatterWorkflowId = WorkflowId.randomId()
  val abortedWorkflowId = WorkflowId.randomId()
}

class MockWorkflowManagerActor extends Actor {
  def receive = {
    case AbortWorkflowCommand(id) =>
      val message = id match {
        case MockWorkflowManagerActor.runningWorkflowId => WorkflowManagerAbortSuccess(id)
        case MockWorkflowManagerActor.abortedWorkflowId =>
          WorkflowManagerAbortFailure(id, new IllegalStateException(s"Workflow ID '$id' is in terminal state 'Aborted' and cannot be aborted."))
      }
      sender ! message
  }
}

class CromwellApiServiceSpec extends FlatSpec with CromwellApiService with ScalatestRouteTest with Matchers
  with ScalaFutures with Mockito {
  import spray.httpx.SprayJsonSupport._

  // BUG: Must be called once to statically initialize the backends, otherwise this Spec won't run if run alone.
  new CromwellSystem {}

  import akka.testkit._

  import scala.concurrent.duration._

  // The submit route takes a bit longer than the default 1 s while things initialize, when this spec is run by itself
  implicit val defaultTimeout = RouteTestTimeout(30.seconds.dilated)

  override def actorRefFactory = system
  override val serviceRegistryActor = CromwellTestkitSpec.ServiceRegistryActorInstance

  override val workflowManagerActor = actorRefFactory.actorOf(Props(new MockWorkflowManagerActor() with WorkflowDescriptorBuilder {
    override implicit  val actorSystem = context.system
  }))

  override val workflowStoreActor = actorRefFactory.actorOf(Props(new MockWorkflowStoreActor()))

  val version = "v1"

  def publishMetadata(events: Seq[MetadataEvent]): Unit = {
    val timeout: Timeout = 5.seconds.dilated

    import akka.pattern.ask
    val putResult = serviceRegistryActor.ask(PutMetadataAction(events))(timeout)
    putResult.futureValue(PatienceConfiguration.Timeout(timeout.duration)) shouldBe a[MetadataPutAcknowledgement]
  }

  def forceSummary(): Unit = {
    val timeout: Timeout = 5.seconds.dilated
    val summaryResult = serviceRegistryActor.ask(RefreshSummary)(timeout)

    val askResult = summaryResult.futureValue(PatienceConfiguration.Timeout(timeout.duration))
    askResult match {
      case MetadataSummarySuccess =>
      case _ =>
        fail()
    }

  }

  behavior of "REST API /status endpoint"

  it should "return 500 errors as Json" in {
    val apiActor = TestActorRef(new CromwellServerActor(ConfigFactory.empty()))
    val probe = TestProbe()
    probe.send(apiActor, Timedout(mock[HttpRequest]))
    probe.expectMsgPF(defaultTimeout.duration) {
      case response: HttpResponse =>
        response.entity.toOption shouldBe defined
        response.entity.toOption.get.contentType.toString() shouldBe ContentTypes.`application/json`.mediaType.value.toString
    }

    system.stop(apiActor)
    system.stop(probe.ref)
  }

  it should "return 404 for get of unknown workflow" in {
    val workflowId = WorkflowId.randomId()

    Get(s"/workflows/$version/$workflowId/status") ~>
      statusRoute ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
      }
  }

   it should "return 400 for get of a malformed workflow id's status" in {
    Get(s"/workflows/$version/foobar/status") ~>
      statusRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        assertResult(
          """{
             |  "status": "fail",
             |  "message": "Invalid workflow ID: 'foobar'."
             |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
  }

  private def publishStatusAndSubmission(workflowId: WorkflowId, state: WorkflowState): Unit = {
    val events = Seq(
      MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.SubmissionTime), MetadataValue(OffsetDateTime.now())),
      MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.Status), MetadataValue(state))
    )
    publishMetadata(events)
  }

  it should "return 200 for get of a known workflow id" in {
    val workflowId = WorkflowId.randomId()
    publishStatusAndSubmission(workflowId, WorkflowSubmitted)

    Get(s"/workflows/$version/$workflowId/status") ~>
      statusRoute ~>
      check {
        status should be(StatusCodes.OK)
        val result = responseAs[JsObject]
        result.fields(WorkflowMetadataKeys.Status) should be(JsString("Submitted"))
      }
  }

  behavior of "REST API /abort endpoint"

  it should "return 404 for abort of unknown workflow" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.unknownId}/abort") ~>
      abortRoute ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
        assertResult(
          s"""{
              |  "status": "fail",
              |  "message": "Unrecognized workflow ID: ${MockWorkflowManagerActor.unknownId}"
              |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
  }

  it should "return 400 for abort of a malformed workflow id" in {
    Post(s"/workflows/$version/foobar/abort") ~>
      abortRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        assertResult(
          """{
              |  "status": "fail",
              |  "message": "Invalid workflow ID: 'foobar'."
              |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
  }

  it should "return 403 for abort of a workflow in a terminal state" in {
    publishStatusAndSubmission(MockWorkflowManagerActor.abortedWorkflowId, WorkflowAborted)

    Post(s"/workflows/$version/${MockWorkflowManagerActor.abortedWorkflowId}/abort") ~>
    abortRoute ~>
    check {
      assertResult(StatusCodes.Forbidden) {
        status
      }
      assertResult(
        s"""{
            |  "status": "error",
            |  "message": "Workflow ID '${MockWorkflowManagerActor.abortedWorkflowId}' is in terminal state 'Aborted' and cannot be aborted."
            |}""".stripMargin
      ) {
        responseAs[String]
      }
    }
  }

  it should "return 200 for abort of a known workflow id" in {
    publishStatusAndSubmission(MockWorkflowManagerActor.runningWorkflowId, WorkflowRunning)

    Post(s"/workflows/$version/${MockWorkflowManagerActor.runningWorkflowId}/abort") ~>
      abortRoute ~>
      check {
        assertResult(
          s"""{
              |  "id": "${MockWorkflowManagerActor.runningWorkflowId.toString}",
              |  "status": "Aborted"
              |}"""
            .stripMargin) {
          responseAs[String]
        }
        assertResult(StatusCodes.OK) {
          status
        }
      }
  }

  behavior of "REST API submission endpoint"
  it should "return 201 for a successful workflow submission " in {
    Post("/workflows/$version", FormData(Seq("wdlSource" -> HelloWorld.wdlSource(), "workflowInputs" -> HelloWorld.rawInputs.toJson.toString()))) ~>
      submitRoute ~>
      check {
        assertResult(
          s"""{
              |  "id": "${MockWorkflowStoreActor.submittedWorkflowId.toString}",
              |  "status": "Submitted"
              |}""".stripMargin) {
          responseAs[String]
        }
        assertResult(StatusCodes.Created) {
          status
        }
      }
  }

  behavior of "REST API batch submission endpoint"
  it should "return 200 for a successful workflow submission " in {
    val inputs = HelloWorld.rawInputs.toJson

    Post("/workflows/$version/batch",
      FormData(Seq("wdlSource" -> HelloWorld.wdlSource(), "workflowInputs" -> s"[$inputs, $inputs]"))) ~>
      submitBatchRoute ~>
      check {
        assertResult(
          s"""[{
              |  "id": "${MockWorkflowStoreActor.submittedWorkflowId.toString}",
              |  "status": "Submitted"
              |}, {
              |  "id": "${MockWorkflowStoreActor.submittedWorkflowId.toString}",
              |  "status": "Submitted"
              |}]""".stripMargin) {
          responseAs[String]
        }
        assertResult(StatusCodes.OK) {
          status
        }
      }
  }

  // TODO: Test tha batch submission returns expected workflow ids in order
  // TODO: Also (assuming we still validate on submit) test a batch of mixed inputs that return submitted and failed

  behavior of "REST API /outputs endpoint"

  it should "return 200 with GET of outputs on successful execution of workflow" in {
    val workflowId = WorkflowId.randomId()
    val events = Seq(
      MetadataEvent(MetadataKey(workflowId, None, s"${WorkflowMetadataKeys.Outputs}:myfirst"), MetadataValue("myOutput"))
    )

    publishMetadata(events)

    Get(s"/workflows/$version/$workflowId/outputs") ~>
      workflowOutputsRoute ~>
      check {
        status should be(StatusCodes.OK)
        val result = responseAs[JsObject]
        result.fields.keys should contain allOf(WorkflowMetadataKeys.Id, WorkflowMetadataKeys.Outputs)
      }
  }

  it should "return 404 with outputs on unknown workflow" in {
    Get(s"/workflows/$version/${MockWorkflowManagerActor.unknownId}/outputs") ~>
    workflowOutputsRoute ~>
    check {
      assertResult(StatusCodes.NotFound) {
        status
      }
    }
  }

  it should "return 405 with POST of outputs on successful execution of workflow" in {
    Post(s"/workflows/$version/${MockWorkflowStoreActor.submittedWorkflowId.toString}/outputs") ~>
      sealRoute(workflowOutputsRoute) ~>
      check {
        assertResult(StatusCodes.MethodNotAllowed) {
          status
        }
      }
  }

  behavior of "REST API /logs endpoint"

  it should "return 200 with paths to stdout/stderr/backend log" in {

    val workflowId = WorkflowId.randomId()
    val jobKey = Option(MetadataJobKey("mycall", None, 1))
    val events = Seq(
      MetadataEvent(MetadataKey(workflowId, jobKey, CallMetadataKeys.Stdout), MetadataValue("stdout.txt")),
      MetadataEvent(MetadataKey(workflowId, jobKey, CallMetadataKeys.Stderr), MetadataValue("stderr.txt")),
      MetadataEvent(MetadataKey(workflowId, jobKey, s"${CallMetadataKeys.BackendLogsPrefix}:log"), MetadataValue("backend.log"))
    )

    publishMetadata(events)

    Get(s"/workflows/$version/$workflowId/logs") ~>
      workflowLogsRoute ~>
      check {
        status should be(StatusCodes.OK)
        val result = responseAs[JsObject]

        val call = result.fields("calls").convertTo[JsObject].fields("mycall").convertTo[Seq[JsObject]].head
        call.fields("stdout") should be(JsString("stdout.txt"))
        call.fields("stderr") should be(JsString("stderr.txt"))
        call.fields("stdout") should be(JsString("stdout.txt"))
        call.fields("backendLogs").convertTo[JsObject].fields("log") should be (JsString("backend.log"))
      }
  }

  it should "return 404 with logs on unknown workflow" in {
    Get(s"/workflows/$version/${MockWorkflowManagerActor.unknownId}/logs") ~>
      workflowLogsRoute ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
      }
  }


  behavior of "REST API /metadata endpoint"

  it should "return with full metadata from the metadata route" in {
    val workflowId = WorkflowId.randomId()
    val events = Seq(
      MetadataEvent(MetadataKey(workflowId, None, "testKey1"), MetadataValue("myValue1")),
      MetadataEvent(MetadataKey(workflowId, None, "testKey2"), MetadataValue("myValue2"))
    )

    publishMetadata(events)

    Get(s"/workflows/$version/$workflowId/metadata") ~>
      metadataRoute ~>
      check {
        status should be(StatusCodes.OK)
        val result = responseAs[JsObject]
        result.fields.keys should contain allOf("testKey1", "testKey2")
        result.fields.keys shouldNot contain("testKey3")
        result.fields("testKey1") should be(JsString("myValue1"))
        result.fields("testKey2") should be(JsString("myValue2"))
      }
  }

  it should "return with included metadata from the metadata route" in {
    val workflowId = WorkflowId.randomId()

    val events = Seq(
      MetadataEvent(MetadataKey(workflowId, None, "testKey1a"), MetadataValue("myValue1a")),
      MetadataEvent(MetadataKey(workflowId, None, "testKey1b"), MetadataValue("myValue1b")),
      MetadataEvent(MetadataKey(workflowId, None, "testKey2a"), MetadataValue("myValue2a")),
      MetadataEvent(MetadataKey(workflowId, None, "testKey2b"), MetadataValue("myValue2b"))
    )

    publishMetadata(events)

    Get(s"/workflows/$version/$workflowId/metadata?includeKey=testKey1&includeKey=testKey2a") ~>
      metadataRoute ~>
      check {
        status should be(StatusCodes.OK)
        val result = responseAs[JsObject]
        result.fields.keys should contain allOf("testKey1a", "testKey1b", "testKey2a")
        result.fields.keys should contain noneOf("testKey2b", "testKey3")
        result.fields("testKey1a") should be(JsString("myValue1a"))
        result.fields("testKey1b") should be(JsString("myValue1b"))
        result.fields("testKey2a") should be(JsString("myValue2a"))
      }
  }

  it should "return with excluded metadata from the metadata route" in {
    val workflowId = WorkflowId.randomId()

    val events = Seq(
      MetadataEvent(MetadataKey(workflowId, None, "testKey1a"), MetadataValue("myValue1a")),
      MetadataEvent(MetadataKey(workflowId, None, "testKey1b"), MetadataValue("myValue1b")),
      MetadataEvent(MetadataKey(workflowId, None, "testKey2a"), MetadataValue("myValue2a")),
      MetadataEvent(MetadataKey(workflowId, None, "testKey2b"), MetadataValue("myValue2b"))
    )

    publishMetadata(events)

    Get(s"/workflows/$version/$workflowId/metadata?excludeKey=testKey2b&excludeKey=testKey3") ~>
      metadataRoute ~>
      check {
        status should be(StatusCodes.OK)
        val result = responseAs[JsObject]
        result.fields.keys should contain allOf("testKey1a", "testKey1b", "testKey2a")
        result.fields.keys should contain noneOf("testKey2b", "testKey3")
        result.fields("testKey1a") should be(JsString("myValue1a"))
        result.fields("testKey1b") should be(JsString("myValue1b"))
        result.fields("testKey2a") should be(JsString("myValue2a"))
      }
  }

  it should "return an error when included and excluded metadata requested from the metadata route" in {
    val workflowId = WorkflowId.randomId()
    publishStatusAndSubmission(workflowId, WorkflowSucceeded)

    Get(s"/workflows/$version/$workflowId/metadata?includeKey=testKey1&excludeKey=testKey2") ~>
      metadataRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        assertResult(
          s"""{
              |  "status": "fail",
              |  "message": "includeKey and excludeKey may not be specified together"
              |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
  }

  behavior of "REST API /timing endpoint"

  it should "return 200 with an HTML document for the timings route" in {
    publishStatusAndSubmission(MockWorkflowStoreActor.submittedWorkflowId, WorkflowSucceeded)

    Get(s"/workflows/$version/${MockWorkflowStoreActor.submittedWorkflowId}/timing") ~>
      timingRoute ~>
      check {
        assertResult(StatusCodes.OK) { status }
        assertResult("<html>") {
          responseAs[String].substring(0, 6)
        }
      }
  }

  behavior of "REST API /query GET endpoint"

  it should "return 400 for a bad query" in {
    Get(s"/workflows/$version/query?BadKey=foo") ~>
      queryRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        assertResult(
          s"""{
              |  "status": "fail",
              |  "message": "Unrecognized query keys: BadKey"
              |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
  }

  it should "return good results for a good query" in {
    val workflowId = WorkflowId.randomId()
    val runningId = WorkflowId.randomId()
    val events = Seq(
      MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.Status), MetadataValue("Succeeded")),
      MetadataEvent(MetadataKey(runningId, None, WorkflowMetadataKeys.Status), MetadataValue("Running"))
    )

    publishMetadata(events)
    forceSummary()

    Get(s"/workflows/$version/query?status=Succeeded&id=$workflowId") ~>
      queryRoute ~>
      check {
        status should be(StatusCodes.OK)
        val results = responseAs[JsObject].fields("results").convertTo[Seq[JsObject]]

        results.head.fields("id") should be(JsString(workflowId.toString))
        results.head.fields("status") should be(JsString("Succeeded"))
      }
  }

  it should "return link headers for pagination when page and pagesize are set for a good query" in {
    val workflowId = WorkflowId.randomId()
    val events = Seq(
      MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.Status), MetadataValue("Succeeded"))
    )

    publishMetadata(events)
    forceSummary()

    Get(s"/workflows/$version/query?status=Succeeded&id=$workflowId&page=1&pagesize=5") ~>
      queryRoute ~>
      check {
        status should be(StatusCodes.OK)
        val results = responseAs[JsObject].fields("results").convertTo[Seq[JsObject]]

        results.head.fields("id") should be(JsString(workflowId.toString))
        results.head.fields("status") should be(JsString("Succeeded"))

        assertResult(4) {
          headers count { header => header.is("link") }
        }
      }
  }

  behavior of "REST API /query POST endpoint"

  it should "return 400 for a bad query map body" in {
    Post(s"/workflows/$version/query", HttpEntity(ContentTypes.`application/json`, """[{"BadKey":"foo"}]""")) ~>
      queryPostRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        assertResult(
          s"""{
              |  "status": "fail",
              |  "message": "Unrecognized query keys: BadKey"
              |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
  }

  it should "return good results for a good query map body" in {
    Post(s"/workflows/$version/query", HttpEntity(ContentTypes.`application/json`, """[{"status":"Succeeded"}]""")) ~>
      queryPostRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult(true) {
          body.asString.contains("\"status\": \"Succeeded\"")
        }
      }
  }

  it should "return link headers for pagination when page and pagesize are set for a good query map body" in {
    Post(s"/workflows/$version/query", HttpEntity(ContentTypes.`application/json`, """[{"status":"Succeeded"},  {"page": "1"}, {"pagesize": "5"}]""")) ~>
      queryPostRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult(true) {
          body.asString.contains("\"status\": \"Succeeded\"")
          (headers count { header => header.is("link") }) == 4
        }
      }
  }


  it should "return good results for a multiple query map body" in {
    Post(s"/workflows/$version/query", HttpEntity(ContentTypes.`application/json`,
      """[{"status":"Succeeded"}, {"status":"Failed"}]""")) ~>
      queryPostRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult(true) {
          body.asString.contains("\"status\": \"Succeeded\"")
        }
      }
  }

  it should "return 400 bad request for a bad query format body" in {
    Post(s"/workflows/$version/query", HttpEntity(ContentTypes.`application/json`, """[{"status":["Succeeded"]}]""")) ~>
      sealRoute(queryPostRoute) ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
      }
  }

  behavior of "REST API /call-caching endpoint"

  ignore should "disallow call caching for a call" taggedAs PostMVP in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching/w.good_call?allow=false") ~>
    callCachingRoute ~>
    check {
      assertResult(StatusCodes.OK) { status }
    }
  }

  ignore should "reject missing 'allow' when disabling call caching for a call" taggedAs PostMVP in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching/w.good_call") ~>
    callCachingRoute ~>
    check {
      assertResult(StatusCodes.BadRequest) { status }
      assertResult(true) {
        responseAs[String].contains("must specify 'allow' exactly once")
      }
    }
  }

  ignore should "reject bogus calls" taggedAs PostMVP in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching/bogus?allow=true") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) { status }
        assertResult(true) {
          responseAs[String].contains("Invalid call")
        }
      }
  }

  ignore should "reject invalid parameter keys when enabling call caching for a call" taggedAs PostMVP in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching/w.good_call?allow=true&bogusKey=foo") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) { status }
        assertResult(true) {
          responseAs[String].contains("Unrecognized parameters: ")
        }
      }
  }

  ignore should "reject bogus workflows when enabling call caching for a call" taggedAs PostMVP in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.unknownId}/call-caching/w.good_call?allow=true") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) { status }
        assertResult(true) {
          responseAs[String].contains("Unknown workflow")
        }
      }
  }

  ignore should "disallow call caching for a workflow" taggedAs PostMVP in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching?allow=false") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.OK) { status }
      }
  }

  ignore should "reject missing 'allow' when disabling call caching for a workflow" taggedAs PostMVP in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) { status }
        assertResult(true) {
          responseAs[String].contains("must specify 'allow' exactly once")
        }
      }
  }

  ignore should "reject invalid parameter keys when enabling call caching for a workflow" taggedAs PostMVP in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching?allow=true&bogusKey=foo") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) { status }
        assertResult(true) {
          responseAs[String].contains("Unrecognized parameters: ")
        }
      }
  }

  ignore should "reject bogus workflows when enabling call caching for a workflow" taggedAs PostMVP in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.unknownId}/call-caching?allow=true") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) { status }
        assertResult(true) {
          responseAs[String].contains("Unknown workflow")
        }
      }
  }
}

