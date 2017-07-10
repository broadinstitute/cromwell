package cromwell.webservice

import akka.actor.{Actor, ActorSystem, Props}
import cromwell.core.{WorkflowId, WorkflowMetadataKeys, WorkflowSubmitted, WorkflowSucceeded}
import akka.http.scaladsl.coding.{Decoder, Gzip}
import akka.http.scaladsl.server.Route
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import spray.json.DefaultJsonProtocol._
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor.{AbortWorkflow, BatchSubmitWorkflows, SubmitWorkflow}
import cromwell.engine.workflow.workflowstore.WorkflowStoreEngineActor
import cromwell.engine.workflow.workflowstore.WorkflowStoreEngineActor.WorkflowAbortFailed
import cromwell.engine.workflow.workflowstore.WorkflowStoreSubmitActor.{WorkflowSubmittedToStore, WorkflowsBatchSubmittedToStore}
import cromwell.services.metadata.MetadataService._
import akka.http.scaladsl.model._
import akka.http.scaladsl.model.headers.{HttpEncodings, `Accept-Encoding`}
import akka.http.scaladsl.testkit.{RouteTestTimeout, ScalatestRouteTest}
import akka.http.scaladsl.unmarshalling.Unmarshal
import akka.stream.ActorMaterializer
import cromwell.services.metadata._
import cromwell.util.SampleWdl.HelloWorld
import org.scalatest.{AsyncFlatSpec, Matchers}
import spray.json._

import scala.concurrent.duration._

class CromwellApiServiceSpec extends AsyncFlatSpec with ScalatestRouteTest with Matchers {
  import CromwellApiServiceSpec._

  val akkaHttpService = new MockApiService()
  val version = "v1"

  implicit def default = RouteTestTimeout(5.seconds)


    behavior of "REST API /status endpoint"
    it should "return 200 for get of a known workflow id" in {
      val workflowId = CromwellApiServiceSpec.ExistingWorkflowId

      Get(s"/workflows/$version/$workflowId/status") ~>
        akkaHttpService.routes ~>
        check {
            status should be(StatusCodes.OK)
            // Along w/ checking value, ensure it is valid JSON despite the requested content type
            responseAs[JsObject].fields(WorkflowMetadataKeys.Status) should be(JsString("Submitted"))
        }
    }

    it should "return 404 for get of unknown workflow" in {
      val workflowId = CromwellApiServiceSpec.UnrecognizedWorkflowId
      Get(s"/workflows/$version/$workflowId/status") ~>
        akkaHttpService.routes ~>
        check {
          assertResult(StatusCodes.NotFound) {
            status
          }
        }
    }

    it should "return 400 for get of a malformed workflow id's status" in {
      Get(s"/workflows/$version/foobar/status") ~>
        akkaHttpService.routes ~>
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

    behavior of "REST API /abort endpoint"
    it should "return 404 for abort of unknown workflow" in {
      val workflowId = CromwellApiServiceSpec.UnrecognizedWorkflowId

      Post(s"/workflows/$version/$workflowId/abort") ~>
        akkaHttpService.routes ~>
        check {
          assertResult(StatusCodes.NotFound) {
            status
          }
        }
    }

    it should "return 400 for abort of a malformed workflow id" in {
      Post(s"/workflows/$version/foobar/abort") ~>
        akkaHttpService.routes ~>
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
      Post(s"/workflows/$version/${CromwellApiServiceSpec.AbortedWorkflowId}/abort") ~>
      akkaHttpService.routes ~>
      check {
        assertResult(StatusCodes.Forbidden) {
          status
        }
        assertResult(
          s"""{
              |  "status": "error",
              |  "message": "Workflow ID '${CromwellApiServiceSpec.AbortedWorkflowId}' is in terminal state 'Aborted' and cannot be aborted."
              |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
    }
    it should "return 200 for abort of a known workflow id" in {
      Post(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/abort") ~>
        akkaHttpService.routes ~>
        check {
          assertResult(
            s"""{"id":"${CromwellApiServiceSpec.ExistingWorkflowId.toString}","status":"Aborted"}""") {
            responseAs[String]
          }
          assertResult(StatusCodes.OK) {
            status
          }
        }
    }

    behavior of "REST API submission endpoint"
    it should "return 201 for a successful workflow submission " in {
      val workflowSource = Multipart.FormData.BodyPart("workflowSource", HttpEntity(MediaTypes.`application/json`, HelloWorld.workflowSource()))
      val workflowInputs = Multipart.FormData.BodyPart("workflowInputs", HttpEntity(MediaTypes.`application/json`, HelloWorld.rawInputs.toJson.toString()))
      val formData = Multipart.FormData(workflowSource, workflowInputs).toEntity()
      Post(s"/workflows/$version", formData) ~>
        akkaHttpService.routes ~>
        check {
          assertResult(
            s"""{
                |  "id": "${CromwellApiServiceSpec.ExistingWorkflowId.toString}",
                |  "status": "Submitted"
                |}""".stripMargin) {
            responseAs[String].parseJson.prettyPrint
          }
          assertResult(StatusCodes.Created) {
            status
          }
        }
    }

    it should "return 400 for an unrecognized form data request parameter " in {
      val formData = Multipart.FormData(Multipart.FormData.BodyPart("incorrectParameter", HttpEntity(MediaTypes.`application/json`, HelloWorld.workflowSource()))).toEntity()
      Post(s"/workflows/$version", formData) ~>
        akkaHttpService.routes ~>
        check {
          assertResult(
            s"""{
                |  "status": "fail",
                |  "message": "Error(s): Unexpected body part name: incorrectParameter"
                |}""".stripMargin) {
            responseAs[String]
          }
          assertResult(StatusCodes.BadRequest) {
            status
          }
        }
    }

    it should "return 400 for a workflow submission with unsupported workflow option keys" in {
      val options = """
                      |{
                      |  "defaultRuntimeOptions": {
                      |  "cpu":1
                      |  }
                      |}
                      |""".stripMargin

      val workflowSource = Multipart.FormData.BodyPart("workflowSource", HttpEntity(MediaTypes.`application/json`, HelloWorld.workflowSource()))
      val workflowInputs = Multipart.FormData.BodyPart("workflowOptions", HttpEntity(MediaTypes.`application/json`, options))
      val formData = Multipart.FormData(workflowSource, workflowInputs).toEntity()

      Post(s"/workflows/$version", formData) ~>
        akkaHttpService.routes ~>
        check {
          assertResult(StatusCodes.BadRequest) {
            status
          }
        }
    }

    it should "return 400 for a workflow submission with malformed workflow options json" in {
      val options = s"""
                       |{"read_from_cache": "true"
                       |""".stripMargin

      val workflowSource = Multipart.FormData.BodyPart("workflowSource", HttpEntity(MediaTypes.`application/json`, HelloWorld.workflowSource()))
      val workflowInputs = Multipart.FormData.BodyPart("workflowOptions", HttpEntity(MediaTypes.`application/json`, options))
      val formData = Multipart.FormData(workflowSource, workflowInputs).toEntity()

      Post(s"/workflows/$version", formData) ~>
        akkaHttpService.routes ~>
        check {
          assertResult(StatusCodes.BadRequest) {
            status
          }
        }
    }

    behavior of "REST API batch submission endpoint"
    it should "return 200 for a successful workflow submission " in {
      val inputs = HelloWorld.rawInputs.toJson
      val workflowSource = Multipart.FormData.BodyPart("workflowSource", HttpEntity(MediaTypes.`application/json`, HelloWorld.workflowSource()))
      val workflowInputs = Multipart.FormData.BodyPart("workflowInputs", HttpEntity(MediaTypes.`application/json`, s"[$inputs, $inputs]"))
      val formData = Multipart.FormData(workflowSource, workflowInputs).toEntity()

      Post(s"/workflows/$version/batch", formData) ~>
        akkaHttpService.routes ~>
        check {
          assertResult(
            s"""[{
                |  "id": "${CromwellApiServiceSpec.ExistingWorkflowId.toString}",
                |  "status": "Submitted"
                |}, {
                |  "id": "${CromwellApiServiceSpec.ExistingWorkflowId.toString}",
                |  "status": "Submitted"
                |}]""".stripMargin) {
            responseAs[String].parseJson.prettyPrint
          }
          assertResult(StatusCodes.Created) {
            status
          }
        }
    }

    it should "return 400 for an submission with no inputs" in {
      val formData = Multipart.FormData(Multipart.FormData.BodyPart("workflowSource", HttpEntity(MediaTypes.`application/json`, HelloWorld.workflowSource()))).toEntity()

      Post(s"/workflows/$version/batch", formData) ~>
       akkaHttpService.routes ~>
        check {
          assertResult(
            s"""{
                |  "status": "fail",
                |  "message": "Error(s): No inputs were provided"
                |}""".stripMargin) {
            responseAs[String]
          }
          assertResult(StatusCodes.BadRequest) {
            status
          }
        }
    }

    behavior of "REST API /outputs endpoint"
    it should "return 200 with GET of outputs on successful execution of workflow" in {
      Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/outputs") ~>
        akkaHttpService.routes ~>
        check {
          status should be(StatusCodes.OK)
          responseAs[JsObject].fields.keys should contain allOf(WorkflowMetadataKeys.Id, WorkflowMetadataKeys.Outputs)
        }
    }

    it should "return 404 with outputs on unknown workflow" in {
      Get(s"/workflows/$version/${CromwellApiServiceSpec.UnrecognizedWorkflowId}/outputs") ~>
        akkaHttpService.routes ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
      }
    }

    it should "return 405 with POST of outputs on successful execution of workflow" in {
      Post(s"/workflows/$version/${CromwellApiServiceSpec.UnrecognizedWorkflowId}/outputs") ~>
        Route.seal(akkaHttpService.routes) ~>
        check {
          assertResult(StatusCodes.MethodNotAllowed) {
            status
          }
        }
    }

    behavior of "REST API /logs endpoint"
    it should "return 200 with paths to stdout/stderr/backend log" in {
      Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/logs") ~>
        akkaHttpService.routes ~>
        check {
          status should be(StatusCodes.OK)

          val call = responseAs[JsObject].fields("calls").convertTo[JsObject].fields("mycall").convertTo[Seq[JsObject]].head
          call.fields("stdout") should be(JsString("stdout.txt"))
          call.fields("stderr") should be(JsString("stderr.txt"))
          call.fields("stdout") should be(JsString("stdout.txt"))
          call.fields("backendLogs").convertTo[JsObject].fields("log") should be (JsString("backend.log"))
        }
    }

    it should "return 404 with logs on unknown workflow" in {
      Get(s"/workflows/$version/${CromwellApiServiceSpec.UnrecognizedWorkflowId}/logs") ~>
        akkaHttpService.routes ~>
        check {
          assertResult(StatusCodes.NotFound) {
            status
          }
        }
    }

    behavior of "REST API /metadata endpoint"
    it should "return with full metadata from the metadata route" in {
      Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/metadata") ~>
        akkaHttpService.routes ~>
        check {
          status should be(StatusCodes.OK)
          val result = responseAs[JsObject]
          result.fields.keys should contain allOf("testKey1", "testKey2")
          result.fields.keys shouldNot contain("testKey3")
          result.fields("testKey1") should be(JsString("myValue1"))
          result.fields("testKey2") should be(JsString("myValue2"))
        }
    }

    it should "return with gzip encoding when requested" in {
      Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/metadata").addHeader(`Accept-Encoding`(HttpEncodings.gzip)) ~>
        akkaHttpService.routes ~>
        check {
          response.headers.find(_.name == "Content-Encoding").get.value should be("gzip")
        }
    }

    it should "not return with gzip encoding when not requested" in {
      Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/metadata") ~>
        akkaHttpService.routes ~>
        check {
          response.headers.find(_.name == "Content-Encoding") shouldBe None
        }
    }

    it should "return with included metadata from the metadata route" in {
      Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/metadata?includeKey=testKey1&includeKey=testKey2a") ~>
       akkaHttpService.routes ~>
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
     Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/metadata?excludeKey=testKey2b&excludeKey=testKey3") ~>
       akkaHttpService.routes ~>
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
      Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/metadata?includeKey=testKey1&excludeKey=testKey2") ~>
        akkaHttpService.routes ~>
        check {
          assertResult(StatusCodes.BadRequest) {
            status
          }

          val decoder: Decoder = Gzip
          Unmarshal(decoder.decodeMessage(response)).to[String] map { r =>
            assertResult(
              s"""{
                  |  "status": "fail",
                  |  "message": "includeKey and excludeKey may not be specified together"
                  |}""".stripMargin
            ) { r }
          }
        }
    }

    behavior of "REST API /timing endpoint"
    it should "return 200 with an HTML document for the timings route" in {
      Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/timing") ~>
        akkaHttpService.routes ~>
        check {
          assertResult(StatusCodes.OK) { status }
          assertResult("<html>") {
            responseAs[String].substring(0, 6)
          }
        }
    }

    behavior of "REST API /query GET endpoint"
    it should "return good results for a good query" in {
      Get(s"/workflows/$version/query?status=Succeeded&id=${CromwellApiServiceSpec.ExistingWorkflowId}") ~>
        akkaHttpService.routes ~>
        check {
          status should be(StatusCodes.OK)
          contentType should be(ContentTypes.`application/json`)
          val results = responseAs[JsObject].fields("results").convertTo[Seq[JsObject]]
          results.head.fields("id") should be(JsString(CromwellApiServiceSpec.ExistingWorkflowId.toString))
          results.head.fields("status") should be(JsString("Succeeded"))
        }
    }

    behavior of "REST API /query POST endpoint"
    it should "return good results for a good query map body" in {
      Post(s"/workflows/$version/query", HttpEntity(ContentTypes.`application/json`, """[{"status":"Succeeded"}]""")) ~>
        akkaHttpService.routes ~>
        check {
          assertResult(StatusCodes.OK) {
            status
          }
          assertResult(true) {
            entityAs[String].contains("\"status\":\"Succeeded\"")
          }
        }
    }

    behavior of "REST API /labels PATCH endpoint"
    it should "return successful status response when assigning valid labels to an existing workflow ID" in {

      val validLabelsJson =
        """
          |{
          |  "label-key-1":"label-value-1",
          |  "label-key-2":"label-value-2"
          |}
        """.stripMargin

      val workflowId = CromwellApiServiceSpec.ExistingWorkflowId

      Patch(s"/workflows/$version/$workflowId/labels", HttpEntity(ContentTypes.`application/json`, validLabelsJson)) ~>
        akkaHttpService.routes ~>
        check {
          status shouldBe StatusCodes.OK
          val actualResult = responseAs[JsObject]
          val expectedResults =
            s"""
              |{
              |  "id": "$workflowId",
              |  "labels": {
              |    "label-key-1":"label-value-1",
              |    "label-key-2":"label-value-2"
              |  }
              |}
            """.stripMargin.parseJson

          actualResult shouldBe expectedResults
        }
    }
}

object CromwellApiServiceSpec {
  val ExistingWorkflowId = WorkflowId.fromString("c4c6339c-8cc9-47fb-acc5-b5cb8d2809f5")
  val AbortedWorkflowId = WorkflowId.fromString("0574111c-c7d3-4145-8190-7a7ed8e8324a")
  val UnrecognizedWorkflowId = WorkflowId.fromString("2bdd06cc-e794-46c8-a897-4c86cedb6a06")
  val RecognizedWorkflowIds = Set(ExistingWorkflowId, AbortedWorkflowId)

  class MockApiService()(implicit val system: ActorSystem) extends CromwellApiService {
    override def actorRefFactory = system

    override val materializer = ActorMaterializer()
    override val ec = system.dispatcher
    override val workflowStoreActor = actorRefFactory.actorOf(Props(new MockWorkflowStoreActor()))
    override val serviceRegistryActor = actorRefFactory.actorOf(Props(new MockServiceRegistryActor()))
    override val workflowManagerActor = actorRefFactory.actorOf(Props.empty)
  }

  object MockServiceRegistryActor {
    def fullMetadataResponse(workflowId: WorkflowId) = {
      List(MetadataEvent(MetadataKey(workflowId, None, "testKey1"), MetadataValue("myValue1", MetadataString)),
        MetadataEvent(MetadataKey(workflowId, None, "testKey2"), MetadataValue("myValue2", MetadataString)))
    }
    def filteredMetadataResponse(workflowId: WorkflowId) = {
      List(MetadataEvent(MetadataKey(workflowId, None, "testKey1a"), MetadataValue("myValue1a", MetadataString)),
        MetadataEvent(MetadataKey(workflowId, None, "testKey1b"), MetadataValue("myValue1b", MetadataString)),
        MetadataEvent(MetadataKey(workflowId, None, "testKey2a"), MetadataValue("myValue2a", MetadataString)))
    }

    def metadataQuery(workflowId: WorkflowId) = MetadataQuery(workflowId, None, None, None, None, false)

    def logsEvents(id: WorkflowId) = {
      val stdout = MetadataEvent(MetadataKey(id, Some(MetadataJobKey("mycall", None, 1)), CallMetadataKeys.Stdout), MetadataValue("stdout.txt", MetadataString))
      val stderr = MetadataEvent(MetadataKey(id, Some(MetadataJobKey("mycall", None, 1)), CallMetadataKeys.Stderr), MetadataValue("stderr.txt", MetadataString))
      val backend = MetadataEvent(MetadataKey(id, Some(MetadataJobKey("mycall", None, 1)), s"${CallMetadataKeys.BackendLogsPrefix}:log"), MetadataValue("backend.log", MetadataString))
      Vector(stdout, stderr, backend)
    }
  }

  class MockServiceRegistryActor extends Actor {
    import MockServiceRegistryActor._
    override def receive = {
      case WorkflowQuery(_) =>
        val response = WorkflowQuerySuccess(WorkflowQueryResponse(List(WorkflowQueryResult(ExistingWorkflowId.toString,
          None, Some(WorkflowSucceeded.toString), None, None))), None)
        sender ! response
      case ValidateWorkflowId(id) =>
        if (RecognizedWorkflowIds.contains(id)) sender ! MetadataService.RecognizedWorkflowId
        else sender ! MetadataService.UnrecognizedWorkflowId
      case GetStatus(id) => sender ! StatusLookupResponse(id, WorkflowSubmitted)
      case WorkflowOutputs(id) =>
        val event = Vector(MetadataEvent(MetadataKey(id, None, "outputs:test.hello.salutation"), MetadataValue("Hello foo!", MetadataString)))
        sender ! WorkflowOutputsResponse(id, event)
      case GetLogs(id) => sender ! LogsResponse(id, logsEvents(id))
      case GetSingleWorkflowMetadataAction(id, None, None, _) => sender ! MetadataLookupResponse(metadataQuery(id), fullMetadataResponse(id))
      case GetSingleWorkflowMetadataAction(id, Some(_), None, _) => sender ! MetadataLookupResponse(metadataQuery(id), filteredMetadataResponse(id))
      case GetSingleWorkflowMetadataAction(id, None, Some(_), _) => sender ! MetadataLookupResponse(metadataQuery(id), filteredMetadataResponse(id))
      case PutMetadataActionAndRespond(events, _) =>
        events.head.key.workflowId match {
          case CromwellApiServiceSpec.ExistingWorkflowId => sender ! MetadataWriteSuccess(events)
          case CromwellApiServiceSpec.AbortedWorkflowId => sender ! MetadataWriteFailure(new Exception("mock exception of db failure"), events)
          case WorkflowId(_) => throw new Exception("Something untoward happened, this situation is not believed to be possible at this time")
        }
    }
  }

  class MockWorkflowStoreActor extends Actor {
    override def receive = {
      case SubmitWorkflow(_) => sender ! WorkflowSubmittedToStore(ExistingWorkflowId)
      case BatchSubmitWorkflows(sources) =>
        val response = WorkflowsBatchSubmittedToStore(sources map { _ => ExistingWorkflowId })
        sender ! response
      case AbortWorkflow(id, manager @ _) =>
        val message = id match {
          case ExistingWorkflowId => WorkflowStoreEngineActor.WorkflowAborted(id)
          case AbortedWorkflowId =>
            WorkflowAbortFailed(id, new IllegalStateException(s"Workflow ID '$id' is in terminal state 'Aborted' and cannot be aborted."))
          case WorkflowId(_) => throw new Exception("Something untoward happened")
        }
        sender ! message
    }
  }
}
