package cromwell.webservice

import akka.actor.{Actor, ActorSystem, Props}
import cromwell.core.{WorkflowId, WorkflowMetadataKeys, WorkflowSubmitted, WorkflowSucceeded}

import scala.util.{Failure, Success, Try}
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor.{AbortWorkflow, BatchSubmitWorkflows, SubmitWorkflow}
import cromwell.engine.workflow.workflowstore.WorkflowStoreEngineActor
import cromwell.engine.workflow.workflowstore.WorkflowStoreEngineActor.WorkflowAbortFailed
import cromwell.services.metadata.MetadataService._
import org.scalatest.{FlatSpec, Matchers}
import spray.http._
import spray.json.DefaultJsonProtocol._
import spray.json.{JsString, _}
import spray.routing._
import cromwell.engine.workflow.workflowstore.WorkflowStoreSubmitActor.{WorkflowSubmittedToStore, WorkflowsBatchSubmittedToStore}
import cromwell.util.SampleWdl.HelloWorld
import spray.httpx.ResponseTransformation
import spray.httpx.SprayJsonSupport._
import spray.httpx.encoding.Gzip
import spray.testkit.ScalatestRouteTest
import spray.routing.Directives._

class CromwellApiServiceSpec extends FlatSpec with ScalatestRouteTest with Matchers with ResponseTransformation {
  import CromwellApiServiceSpec._

  val cromwellApiService = new MockApiService()
  val version = "v1"

  behavior of "REST API /status endpoint"
    it should "return 200 for get of a known workflow id" in {
      val workflowId = MockApiService.ExistingWorkflowId
      Get(s"/workflows/$version/$workflowId/status") ~>
        cromwellApiService.statusRoute ~>
        check {
            status should be(StatusCodes.OK)
            val result = responseAs[JsObject]
            result.fields(WorkflowMetadataKeys.Status) should be(JsString("Submitted"))
        }
    }

    it should "return 404 for get of unknown workflow" in {
      val workflowId = MockApiService.UnrecognizedWorkflowId

      Get(s"/workflows/$version/$workflowId/status") ~>
        cromwellApiService.statusRoute ~>
        check {
          assertResult(StatusCodes.NotFound) {
            status
          }
        }
    }

    it should "return 400 for get of a malformed workflow id's status" in {
      Get(s"/workflows/$version/foobar/status") ~>
        cromwellApiService.statusRoute ~>
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
      val workflowId = MockApiService.UnrecognizedWorkflowId

      Post(s"/workflows/$version/$workflowId/abort") ~>
        cromwellApiService.abortRoute ~>
        check {
          assertResult(StatusCodes.NotFound) {
            status
          }
        }
    }

    it should "return 400 for abort of a malformed workflow id" in {
      Post(s"/workflows/$version/foobar/abort") ~>
        cromwellApiService.abortRoute ~>
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
      Post(s"/workflows/$version/${MockApiService.AbortedWorkflowId}/abort") ~>
      cromwellApiService.abortRoute ~>
      check {
        assertResult(StatusCodes.Forbidden) {
          status
        }
        assertResult(
          s"""{
              |  "status": "error",
              |  "message": "Workflow ID '${MockApiService.AbortedWorkflowId}' is in terminal state 'Aborted' and cannot be aborted."
              |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
    }

    it should "return 200 for abort of a known workflow id" in {
      Post(s"/workflows/$version/${MockApiService.ExistingWorkflowId}/abort") ~>
        cromwellApiService.abortRoute ~>
        check {
          assertResult(
            s"""{
                |  "id": "${MockApiService.ExistingWorkflowId.toString}",
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
    val bodyParts: Map[String, BodyPart] = Map("wdlSource" -> BodyPart(HelloWorld.wdlSource()), "workflowInputs" -> BodyPart(HelloWorld.rawInputs.toJson.toString()))
    Post(s"/workflows/$version", MultipartFormData(bodyParts)) ~>
      cromwellApiService.submitRoute ~>
      check {
        assertResult(
          s"""{
              |  "id": "${MockApiService.ExistingWorkflowId.toString}",
              |  "status": "Submitted"
              |}""".stripMargin) {
          responseAs[String]
        }
        assertResult(StatusCodes.Created) {
          status
        }
      }
  }

  it should "return 400 for an unrecognized form data request parameter " in {
    val bodyParts: Map[String, BodyPart] = Map("incorrectParameter" -> BodyPart(HelloWorld.wdlSource()))
    Post(s"/workflows/$version", MultipartFormData(bodyParts)) ~>
      cromwellApiService.submitRoute ~>
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

    val bodyParts = Map("wdlSource" -> BodyPart(HelloWorld.wdlSource()), "workflowOptions" -> BodyPart(options))

    Post(s"/workflows/$version", MultipartFormData(bodyParts)) ~>
      cromwellApiService.submitRoute ~>
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

    val bodyParts = Map("wdlSource" -> BodyPart(HelloWorld.wdlSource()), "workflowOptions" -> BodyPart(options))

    Post(s"/workflows/$version", MultipartFormData(bodyParts)) ~>
      cromwellApiService.submitRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
      }
  }

  it should "succesfully merge and override multiple input files" in {

    val input1 = Map("wf.a1" -> "hello", "wf.a2" -> "world").toJson.toString
    val input2 = Map.empty[String, String].toJson.toString
    val overrideInput1 = Map("wf.a2" -> "universe").toJson.toString
    val allInputs = cromwellApiService.mergeMaps(Seq(Option(input1), Option(input2), Option(overrideInput1)))

    check {
      allInputs.fields.keys should contain allOf("wf.a1", "wf.a2")
      allInputs.fields("wf.a2") should be(JsString("universe"))
    }
  }

  behavior of "REST API batch submission endpoint"
  it should "return 200 for a successful workflow submission " in {
    val inputs = HelloWorld.rawInputs.toJson
    val bodyParts = Map("wdlSource" -> BodyPart(HelloWorld.wdlSource()), "workflowInputs" -> BodyPart(s"[$inputs, $inputs]"))

    Post(s"/workflows/$version/batch", MultipartFormData(bodyParts)) ~>
      cromwellApiService.submitBatchRoute ~>
      check {
        assertResult(
          s"""[{
              |  "id": "${MockApiService.ExistingWorkflowId.toString}",
              |  "status": "Submitted"
              |}, {
              |  "id": "${MockApiService.ExistingWorkflowId.toString}",
              |  "status": "Submitted"
              |}]""".stripMargin) {
          responseAs[String]
        }
        assertResult(StatusCodes.OK) {
          status
        }
      }
  }

  it should "return 400 for an submission with no inputs" in {
    val bodyParts = Map("wdlSource" -> BodyPart(HelloWorld.wdlSource()))

    Post(s"/workflows/$version/batch", MultipartFormData(bodyParts)) ~>
      cromwellApiService.submitBatchRoute ~>
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
    Get(s"/workflows/$version/${MockApiService.ExistingWorkflowId}/outputs") ~>
      cromwellApiService.workflowOutputsRoute ~>
      check {
        status should be(StatusCodes.OK)
        val result = responseAs[JsObject]
        result.fields.keys should contain allOf(WorkflowMetadataKeys.Id, WorkflowMetadataKeys.Outputs)
      }
  }

  it should "return 404 with outputs on unknown workflow" in {
    Get(s"/workflows/$version/${MockApiService.UnrecognizedWorkflowId}/outputs") ~>
    cromwellApiService.workflowOutputsRoute ~>
    check {
      assertResult(StatusCodes.NotFound) {
        status
      }
    }
  }

  it should "return 405 with POST of outputs on successful execution of workflow" in {
    Post(s"/workflows/$version/${MockApiService.UnrecognizedWorkflowId}/outputs") ~>
      cromwellApiService.sealRoute(cromwellApiService.workflowOutputsRoute) ~>
      check {
        assertResult(StatusCodes.MethodNotAllowed) {
          status
        }
      }
  }

  behavior of "REST API /logs endpoint"
  it should "return 200 with paths to stdout/stderr/backend log" in {
    Get(s"/workflows/$version/${MockApiService.ExistingWorkflowId}/logs") ~>
      cromwellApiService.workflowLogsRoute ~>
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
    Get(s"/workflows/$version/${MockApiService.UnrecognizedWorkflowId}/logs") ~>
      cromwellApiService.workflowLogsRoute ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
      }
  }

  behavior of "REST API /metadata endpoint"
  it should "return with full metadata from the metadata route" in {
    Get(s"/workflows/$version/${MockApiService.ExistingWorkflowId}/metadata") ~>
      mapHttpResponse(decode(Gzip))(cromwellApiService.metadataRoute) ~>
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
    Get(s"/workflows/$version/${MockApiService.ExistingWorkflowId}/metadata?includeKey=testKey1&includeKey=testKey2a") ~>
      mapHttpResponse(decode(Gzip))(cromwellApiService.metadataRoute) ~>
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
     Get(s"/workflows/$version/${MockApiService.ExistingWorkflowId}/metadata?excludeKey=testKey2b&excludeKey=testKey3") ~>
       mapHttpResponse(decode(Gzip))(cromwellApiService.metadataRoute) ~>
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
    Get(s"/workflows/$version/${MockApiService.ExistingWorkflowId}/metadata?includeKey=testKey1&excludeKey=testKey2") ~>
      mapHttpResponse(decode(Gzip))(cromwellApiService.metadataRoute) ~>
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
    Get(s"/workflows/$version/${MockApiService.ExistingWorkflowId}/timing") ~>
      cromwellApiService.timingRoute ~>
      check {
        assertResult(StatusCodes.OK) { status }
        assertResult("<html>") {
          responseAs[String].substring(0, 6)
        }
      }
  }

  behavior of "REST API /query GET endpoint"
  it should "return good results for a good query" in {
    Get(s"/workflows/$version/query?status=Succeeded&id=${MockApiService.ExistingWorkflowId}") ~>
      cromwellApiService.queryRoute ~>
      check {
        status should be(StatusCodes.OK)
        val results = responseAs[JsObject].fields("results").convertTo[Seq[JsObject]]

        results.head.fields("id") should be(JsString(MockApiService.ExistingWorkflowId.toString))
        results.head.fields("status") should be(JsString("Succeeded"))
      }
  }

  behavior of "REST API /query POST endpoint"
  it should "return good results for a good query map body" in {
    Post(s"/workflows/$version/query", HttpEntity(ContentTypes.`application/json`, """[{"status":"Succeeded"}]""")) ~>
      cromwellApiService.queryPostRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult(true) {
          body.asString.contains("\"status\": \"Succeeded\"")
        }
      }
  }
}

object CromwellApiServiceSpec {
  class MockApiService()(implicit val system: ActorSystem) extends CromwellApiService {
    import MockApiService._

    override def actorRefFactory = system
    override val workflowStoreActor = actorRefFactory.actorOf(Props(new MockWorkflowStoreActor()))
    override val serviceRegistryActor = actorRefFactory.actorOf(Props.empty)
    override val workflowManagerActor = actorRefFactory.actorOf(Props.empty)
    override val callCacheReadActor = actorRefFactory.actorOf(Props.empty)

    override def handleMetadataRequest(message: AnyRef): Route = {
      message match {
        case GetStatus(w) => complete(submittedStatusResponse(w))
        case WorkflowOutputs(w) => complete(outputResponse(w))
        case GetLogs(w) => complete(logsResponse(w))
        case GetSingleWorkflowMetadataAction(w, None, None, _) => complete(fullMetadataResponse(w))
        case GetSingleWorkflowMetadataAction(w, Some(i), None, _) => complete(filteredMetadataResponse(w))
        case GetSingleWorkflowMetadataAction(w, None, Some(_), _) => complete(filteredMetadataResponse(w))
        case _ => throw new IllegalArgumentException("Woopsie!")
      }
    }

    override def handleQueryMetadataRequest(parameters: Seq[(String, String)]): Route = {
      complete(QueryResponse)
    }

    override def withRecognizedWorkflowId(possibleWorkflowId: String)(recognizedWorkflowId: WorkflowId => Route): Route = {
      requestContext =>
        Try(WorkflowId.fromString(possibleWorkflowId)) match {
          case Success(workflowId) =>
              if (RecognizedWorkflowIds.contains(workflowId)) recognizedWorkflowId(workflowId)(requestContext)
              else failBadRequest(new NoSuchElementException("A hollow voice says 'fool'"), StatusCodes.NotFound)(requestContext)
          case Failure(e) =>
            failBadRequest(new RuntimeException(s"Invalid workflow ID: '$possibleWorkflowId'."))(requestContext)
        }
    }
  }

  object MockApiService {
    val ExistingWorkflowId = WorkflowId.fromString("c4c6339c-8cc9-47fb-acc5-b5cb8d2809f5")
    val AbortedWorkflowId = WorkflowId.fromString("0574111c-c7d3-4145-8190-7a7ed8e8324a")
    val UnrecognizedWorkflowId = WorkflowId.fromString("2bdd06cc-e794-46c8-a897-4c86cedb6a06")
    val RecognizedWorkflowIds = Set(ExistingWorkflowId, AbortedWorkflowId)

    val QueryResponse = JsObject(Map(
      "results" -> JsArray(JsObject(Map(
        WorkflowMetadataKeys.Id -> JsString(ExistingWorkflowId.toString),
        WorkflowMetadataKeys.Status -> JsString(WorkflowSucceeded.toString)
      )))
    ))

    def fullMetadataResponse(workflowId: WorkflowId) = JsObject(Map(
      "testKey1" -> JsString("myValue1"),
      "testKey2" -> JsString("myValue2")
    ))

    def filteredMetadataResponse(workflowId: WorkflowId) = JsObject(Map(
      "testKey1a" -> JsString("myValue1a"),
      "testKey1b" -> JsString("myValue1b"),
      "testKey2a" -> JsString("myValue2a")
    ))

    def submittedStatusResponse(workflowId: WorkflowId) = JsObject(Map(
      WorkflowMetadataKeys.Status -> JsString(WorkflowSubmitted.toString),
      WorkflowMetadataKeys.Id -> JsString(workflowId.toString)
    ))

    def outputResponse(workflowId: WorkflowId) = JsObject(Map(
      WorkflowMetadataKeys.Id -> JsString(workflowId.toString),
      WorkflowMetadataKeys.Outputs -> JsString("Some random stuff")
    ))

    def logsResponse(workflowId: WorkflowId) = JsObject(Map(
      "calls" -> JsObject(Map(
        "mycall" -> JsArray(JsObject(Map(
          "stdout" -> JsString("stdout.txt"),
          "stderr" -> JsString("stderr.txt"),
          "backendLogs" -> JsObject(Map(
            "log" -> JsString("backend.log")
          ))
        )))
      ))
    ))
  }

  class MockWorkflowStoreActor extends Actor {
    override def receive = {
      case SubmitWorkflow(source) => sender ! WorkflowSubmittedToStore(MockApiService.ExistingWorkflowId)
      case BatchSubmitWorkflows(sources) =>
        val response = WorkflowsBatchSubmittedToStore(sources map { _ => MockApiService.ExistingWorkflowId })
        sender ! response
      case AbortWorkflow(id, manager) =>
        val message = id match {
          case MockApiService.ExistingWorkflowId =>
            WorkflowStoreEngineActor.WorkflowAborted(id)
          case MockApiService.AbortedWorkflowId   =>
            WorkflowAbortFailed(id, new IllegalStateException(s"Workflow ID '$id' is in terminal state 'Aborted' and cannot be aborted."))
        }

        sender ! message
    }
  }
}
