package cromwell.webservice.routes

import akka.actor.{Actor, ActorLogging, ActorSystem, Props}
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model._
import akka.http.scaladsl.testkit.{RouteTestTimeout, ScalatestRouteTest}
import akka.stream.ActorMaterializer
import common.util.VersionUtil
import cromwell.core._
import cromwell.core.abort.{WorkflowAbortFailureResponse, WorkflowAbortRequestedResponse}
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowNotFoundException
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor.{AbortWorkflowCommand, BatchSubmitWorkflows, SubmitWorkflow, WorkflowOnHoldToSubmittedCommand}
import cromwell.engine.workflow.workflowstore.WorkflowStoreEngineActor.{WorkflowOnHoldToSubmittedFailure, WorkflowOnHoldToSubmittedSuccess}
import cromwell.engine.workflow.workflowstore.WorkflowStoreSubmitActor.{WorkflowSubmittedToStore, WorkflowsBatchSubmittedToStore}
import cromwell.services.healthmonitor.HealthMonitorServiceActor.{GetCurrentStatus, StatusCheckResponse, SubsystemStatus}
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import cromwell.util.SampleWdl.HelloWorld
import cromwell.webservice.EngineStatsActor
import mouse.boolean._
import org.scalatest.{AsyncFlatSpec, Matchers}
import spray.json.DefaultJsonProtocol._
import spray.json._

import scala.concurrent.duration._

class CromwellApiServiceSpec extends AsyncFlatSpec with ScalatestRouteTest with Matchers {
  import CromwellApiServiceSpec._

  val akkaHttpService = new MockApiService()
  val version = "v1"

  implicit def default = RouteTestTimeout(5.seconds)

  "REST ENGINE /stats endpoint" should "return 200 for stats" ignore {
    Get(s"/engine/$version/stats") ~>
      akkaHttpService.engineRoutes ~>
      check {
        status should be(StatusCodes.OK)
        val resp = responseAs[JsObject]
        val workflows = resp.fields("workflows").asInstanceOf[JsNumber].value.toInt
        val jobs = resp.fields("jobs").asInstanceOf[JsNumber].value.toInt
        workflows should be(1)
        jobs should be(23)
      }
  }

  "REST ENGINE /version endpoint" should "return 200 for version" in {
    Get(s"/engine/$version/version") ~>
      akkaHttpService.engineRoutes ~>
      check {
        status should be(StatusCodes.OK)
        val resp = responseAs[JsObject]
        val cromwellVersion = resp.fields("cromwell").asInstanceOf[JsString].value
        cromwellVersion should fullyMatch regex
          s"""(\\d+-([0-9a-f]){7}(-SNAP)?|${VersionUtil.defaultMessage("cromwell-engine")})"""
      }
  }

  "REST ENGINE /status endpoint" should "return 200 for status when all is well" in {
    Get(s"/engine/$version/status") ~>
      akkaHttpService.engineRoutes ~>
        check {
          status should be(StatusCodes.OK)
          val resp = responseAs[JsObject]
          val db = resp.fields("Engine Database").asJsObject
          db.fields("ok").asInstanceOf[JsBoolean].value should be(true)
        }
  }


    behavior of "REST API /abort endpoint"
    it should "return 404 for abort of unknown workflow" in {
      val workflowId = CromwellApiServiceSpec.UnrecognizedWorkflowId

      Post(s"/workflows/$version/$workflowId/abort") ~>
        akkaHttpService.workflowRoutes ~>
        check {
          assertResult(StatusCodes.NotFound) {
            status
          }
        }
    }

    it should "return 400 for abort of a malformed workflow id" in {
      Post(s"/workflows/$version/foobar/abort") ~>
        akkaHttpService.workflowRoutes ~>
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
        akkaHttpService.workflowRoutes ~>
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
        akkaHttpService.workflowRoutes ~>
        check {
          assertResult(
            s"""{"id":"${CromwellApiServiceSpec.ExistingWorkflowId.toString}","status":"Aborting"}""") {
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
        akkaHttpService.workflowRoutes ~>
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
          headers should be(Seq.empty)
        }
    }

  it should "return 201 for a successful workflow submission using workflowUrl" in {
    val workflowUrl = Multipart.FormData.BodyPart("workflowUrl", HttpEntity("https://raw.githubusercontent.com/broadinstitute/cromwell/develop/womtool/src/test/resources/validate/wdl_draft3/valid/callable_imports/my_workflow.wdl"))
    val formData = Multipart.FormData(workflowUrl).toEntity()

    Post(s"/workflows/$version", formData) ~>
      akkaHttpService.workflowRoutes ~>
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
        headers should be(Seq.empty)
      }
  }

  it should "return 400 for a workflow submission using workflowUrl with invalid protocol" in {
    val workflowUrl = Multipart.FormData.BodyPart("workflowUrl", HttpEntity("htpps://raw.githubusercontent.com/broadinstitute/cromwell/develop/womtool/src/test/resources/validate/wdl_draft3/valid/callable_imports/my_workflow.wdl"))
    val formData = Multipart.FormData(workflowUrl).toEntity()

    Post(s"/workflows/$version", formData) ~>
      akkaHttpService.workflowRoutes ~>
      check {
        assertResult(
          s"""{
             |  "status": "fail",
             |  "message": "Error(s): Error while validating workflow url: unknown protocol: htpps"
             |}""".stripMargin) {
          responseAs[String].parseJson.prettyPrint
        }
        assertResult(StatusCodes.BadRequest) {
          status
        }
        headers should be(Seq.empty)
      }
  }

  it should "return 201 for a successful workflow submission with onHold = true" in {
    val workflowSource = Multipart.FormData.BodyPart("workflowSource", HttpEntity(MediaTypes.`application/json`, HelloWorld.workflowSource()))
    val workflowInputs = Multipart.FormData.BodyPart("workflowInputs", HttpEntity(MediaTypes.`application/json`, HelloWorld.rawInputs.toJson.toString()))
    val onHold =  Multipart.FormData.BodyPart("workflowOnHold", HttpEntity("true"))
    val formData = Multipart.FormData(workflowSource, workflowInputs, onHold).toEntity()
    Post(s"/workflows/$version", formData) ~>
      akkaHttpService.workflowRoutes ~>
      check {
        assertResult(
          s"""{
             |  "id": "${CromwellApiServiceSpec.ExistingWorkflowId.toString}",
             |  "status": "On Hold"
             |}""".stripMargin) {
          responseAs[String].parseJson.prettyPrint
        }
        assertResult(StatusCodes.Created) {
          status
        }
        headers should be(Seq.empty)
      }
  }

  it should "return 200 when a workflow is switched from on hold to submitted" in {
    val id = ExistingWorkflowId
    Post(s"/workflows/$version/$id/releaseHold") ~>
      akkaHttpService.workflowRoutes ~>
      check {
        assertResult(
          s"""{
             |  "id": "${CromwellApiServiceSpec.ExistingWorkflowId.toString}",
             |  "status": "Submitted"
             |}""".stripMargin) {
          responseAs[String].parseJson.prettyPrint
        }
        assertResult(StatusCodes.OK) {
          status
        }
        headers should be(Seq.empty)
      }
  }

  it should "return 404 when invalid workflow id is submitted to onHoldToSubmitted API end point" in {
    val id = UnrecognizedWorkflowId
    Post(s"/workflows/$version/$id/releaseHold") ~>
      akkaHttpService.workflowRoutes ~>
      check {
        assertResult(
          s"""{
             |  "status": "fail",
             |  "message": "Unrecognized workflow ID: ${CromwellApiServiceSpec.UnrecognizedWorkflowId.toString}"
             |}""".stripMargin) {
          responseAs[String].parseJson.prettyPrint
        }
        assertResult(StatusCodes.NotFound) {
          status
        }
        headers should be(Seq.empty)
      }
  }

  it should "return 201 with warnings for a successful v1 workflow submission still using wdlSource" in {
    val workflowSource = Multipart.FormData.BodyPart("wdlSource",
      HttpEntity(MediaTypes.`application/json`, HelloWorld.workflowSource()))
    val formData = Multipart.FormData(workflowSource).toEntity()
    Post(s"/workflows/v1", formData) ~>
      akkaHttpService.workflowRoutes ~>
      check {
        status should be(StatusCodes.Created)
        responseAs[String].parseJson.prettyPrint should be(
          s"""|{
              |  "id": "${CromwellApiServiceSpec.ExistingWorkflowId.toString}",
              |  "status": "Submitted"
              |}
              |""".stripMargin.trim
        )
        headers.size should be(1)
        val warningHeader = header("Warning")
        warningHeader shouldNot be(empty)
        warningHeader.get.value should fullyMatch regex
          s"""299 cromwell/(\\d+-([0-9a-f]){7}(-SNAP)?|${VersionUtil.defaultMessage("cromwell-engine")}) """ +
            "\"The 'wdlSource' parameter name has been deprecated in favor of 'workflowSource'. " +
            "Support for 'wdlSource' will be removed from future versions of Cromwell. " +
            "Please switch to using 'workflowSource' in future submissions.\""
      }
  }

    it should "return 400 for an unrecognized form data request parameter " in {
      val formData = Multipart.FormData(Map(
        "incorrectParameter" -> HttpEntity(MediaTypes.`application/json`, HelloWorld.workflowSource()),
        "incorrectParameter2" -> HttpEntity(MediaTypes.`application/json`, HelloWorld.workflowSource())
      )).toEntity

      Post(s"/workflows/$version", formData) ~>
        akkaHttpService.workflowRoutes ~>
        check {
          assertResult(
            s"""{
                |  "status": "fail",
                |  "message": "Error(s): Unexpected body part name: incorrectParameter\\nUnexpected body part name: incorrectParameter2\\nworkflowSource or workflowUrl needs to be supplied"
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
        akkaHttpService.workflowRoutes ~>
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
        akkaHttpService.workflowRoutes ~>
        check {
          assertResult(StatusCodes.BadRequest) {
            status
          }
        }
    }

    it should "return 400 for a workflow submission with invalid workflow custom labels" in {
      val labels = s"""
                      |{"key with more than 255 characters-at vero eos et accusamus et iusto odio dignissimos ducimus qui blanditiis praesentium voluptatum deleniti atque corrupti quos dolores et quas molestias excepturi sint occaecati cupiditate non provident, similique sunt in culpas":"value with more than 255 characters-at vero eos et accusamus et iusto odio dignissimos ducimus qui blanditiis praesentium voluptatum deleniti atque corrupti quos dolores et quas molestias excepturi sint occaecati cupiditate non provident, similique sunt in culpa"}
                      |""".stripMargin

      val workflowSource = Multipart.FormData.BodyPart("workflowSource", HttpEntity(MediaTypes.`application/json`, HelloWorld.workflowSource()))
      val customLabels =  Multipart.FormData.BodyPart("labels", HttpEntity(MediaTypes.`application/json`, labels))
      val onHold =  Multipart.FormData.BodyPart("workflowOnHold", HttpEntity("true"))
      val formData = Multipart.FormData(workflowSource, customLabels, onHold).toEntity()

      Post(s"/workflows/$version", formData) ~>
        akkaHttpService.workflowRoutes ~>
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
        akkaHttpService.workflowRoutes ~>
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
        akkaHttpService.workflowRoutes ~>
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

    behavior of "REST API /timing endpoint"
    it should "return 200 with an HTML document for the timings route" in {
      Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/timing") ~>
        akkaHttpService.workflowRoutes ~>
        check {
          assertResult(StatusCodes.OK) { status }
          assertResult("<html>") {
            responseAs[String].substring(0, 6)
          }
        }
    }

    it should "return 200 with an HTML document that contains static metadata json" in {
      Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/timing") ~>
        akkaHttpService.workflowRoutes ~>
        check {
          assertResult(StatusCodes.OK) { status }
          assert(responseAs[String].contains("var metadataJson = {\"testKey1b\":\"myValue1b\",\"calls\":{},\"testKey1a\":\"myValue1a\",\"id\":\"c4c6339c-8cc9-47fb-acc5-b5cb8d2809f5\",\"testKey2a\":\"myValue2a\"};"))
        }
    }

    it should "return 404 when unrecognized workflow id is submitted" in {
      Get(s"/workflows/$version/${CromwellApiServiceSpec.UnrecognizedWorkflowId}/timing") ~>
        akkaHttpService.workflowRoutes ~>
        check {
          assertResult(StatusCodes.NotFound) { status }
          assertResult(
            s"""{
               |  "status": "fail",
               |  "message": "Unrecognized workflow ID: ${CromwellApiServiceSpec.UnrecognizedWorkflowId.toString}"
               |}""".stripMargin) {
            responseAs[String].parseJson.prettyPrint
          }
        }
    }

    it should "return 400 when invalid workflow id is submitted" in {
      Get(s"/workflows/$version/foo/timing") ~>
        akkaHttpService.workflowRoutes ~>
        check {
          assertResult(StatusCodes.BadRequest) { status }
          assertResult(
            s"""{
               |  "status": "fail",
               |  "message": "Invalid workflow ID: 'foo'."
               |}""".stripMargin) {
            responseAs[String].parseJson.prettyPrint
          }
        }
    }
}

object CromwellApiServiceSpec {
  val ExistingWorkflowId = WorkflowId.fromString("c4c6339c-8cc9-47fb-acc5-b5cb8d2809f5")
  val AbortedWorkflowId = WorkflowId.fromString("0574111c-c7d3-4145-8190-7a7ed8e8324a")
  val UnrecognizedWorkflowId = WorkflowId.fromString("2bdd06cc-e794-46c8-a897-4c86cedb6a06")
  val RecognizedWorkflowIds = Set(ExistingWorkflowId)

  class MockApiService()(implicit val system: ActorSystem) extends CromwellApiService {
    override def actorRefFactory = system

    override val materializer = ActorMaterializer()
    override val ec = system.dispatcher
    override val workflowStoreActor = actorRefFactory.actorOf(Props(new MockWorkflowStoreActor()))
    override val serviceRegistryActor = actorRefFactory.actorOf(Props(new MockServiceRegistryActor()))
    override val workflowManagerActor = actorRefFactory.actorOf(Props(new MockWorkflowManagerActor()))
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

    def metadataQuery(workflowId: WorkflowId) =
      MetadataQuery(workflowId, None, None, None, None, expandSubWorkflows = false)

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
      case WorkflowQuery(parameters) =>
        val labels: Option[Map[String, String]] = {
          parameters.contains(("additionalQueryResultFields", "labels")).option(
            Map("key1" -> "label1", "key2" -> "label2"))
        }
        val parentWorkflowId: Option[String] =  {
          parameters.contains(("additionalQueryResultFields", "parentWorkflowId")).option("pid")
        }
        val response = WorkflowQuerySuccess(WorkflowQueryResponse(List(WorkflowQueryResult(ExistingWorkflowId.toString,
          None, Some(WorkflowSucceeded.toString), None, None, None, labels, parentWorkflowId)), 1), None)
        sender ! response
      case ValidateWorkflowId(id) =>
        if (RecognizedWorkflowIds.contains(id)) sender ! MetadataService.RecognizedWorkflowId
        else sender ! MetadataService.UnrecognizedWorkflowId
      case GetCurrentStatus =>
        sender ! StatusCheckResponse(
          ok = true,
          systems = Map(
            "Engine Database" -> SubsystemStatus(ok = true, messages = None)))
      case GetStatus(id) => sender ! StatusLookupResponse(id, WorkflowSubmitted)
      case GetLabels(id) => sender ! LabelLookupResponse(id, Map("key1" -> "label1", "key2" -> "label2"))
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
      case command: WorkflowOnHoldToSubmittedCommand if command.id == ExistingWorkflowId =>
        sender ! WorkflowOnHoldToSubmittedSuccess(command.id)
      case command: WorkflowOnHoldToSubmittedCommand if command.id == UnrecognizedWorkflowId =>
        sender ! WorkflowOnHoldToSubmittedFailure(command.id, new Exception("Cannot switch to submitted"))
      case SubmitWorkflow(_) => sender ! WorkflowSubmittedToStore(ExistingWorkflowId, WorkflowSubmitted)
      case BatchSubmitWorkflows(sources) =>
        val response = WorkflowsBatchSubmittedToStore(sources map { _ => ExistingWorkflowId }, WorkflowSubmitted)
        sender ! response
      case AbortWorkflowCommand(id) =>
        val message = id match {
          case ExistingWorkflowId => WorkflowAbortRequestedResponse(id)
          case UnrecognizedWorkflowId => WorkflowAbortFailureResponse(id, new WorkflowNotFoundException(s"Couldn't abort $id because no workflow with that ID is in progress"))
          case AbortedWorkflowId =>
            WorkflowAbortFailureResponse(id, new IllegalStateException(s"Workflow ID '$id' is in terminal state 'Aborted' and cannot be aborted."))
          case WorkflowId(_) => throw new Exception("Something untoward happened")
        }
        sender ! message
    }
  }
  
  class MockWorkflowManagerActor extends Actor with ActorLogging {
    override def receive: Receive = {
      case WorkflowManagerActor.EngineStatsCommand =>
        val response = EngineStatsActor.EngineStats(1, 23)
        sender ! response
      case unexpected =>
        val sndr = sender()
        log.error(s"Unexpected message {} from {}", unexpected, sndr)
        sender ! s"Unexpected message received: $unexpected"
    }
  }
}
