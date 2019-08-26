package cromwell.webservice.routes

import akka.actor.{Actor, ActorLogging, ActorSystem, Props}
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model.ContentTypes._
import akka.http.scaladsl.model._
import akka.http.scaladsl.testkit.{RouteTestTimeout, ScalatestRouteTest}
import akka.stream.ActorMaterializer
import com.typesafe.scalalogging.StrictLogging
import common.util.VersionUtil
import cromwell.core._
import cromwell.core.abort.{WorkflowAbortFailureResponse, WorkflowAbortRequestedResponse, WorkflowAbortedResponse}
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowNotFoundException
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor._
import cromwell.engine.workflow.workflowstore.WorkflowStoreEngineActor.{WorkflowOnHoldToSubmittedFailure, WorkflowOnHoldToSubmittedSuccess}
import cromwell.engine.workflow.workflowstore.WorkflowStoreSubmitActor.{WorkflowSubmittedToStore, WorkflowsBatchSubmittedToStore}
import cromwell.services.healthmonitor.ProtoHealthMonitorServiceActor.{GetCurrentStatus, StatusCheckResponse, SubsystemStatus}
import cromwell.services.instrumentation.InstrumentationService.InstrumentationServiceMessage
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import cromwell.services.metadata.impl.builder.MetadataBuilderActor
import cromwell.services.metadata.impl.builder.MetadataBuilderActor.BuiltMetadataResponse
import cromwell.services.womtool.WomtoolServiceMessages.{DescribeFailure, DescribeRequest, DescribeSuccess}
import cromwell.services.womtool.models.WorkflowDescription
import cromwell.util.SampleWdl.HelloWorld
import cromwell.webservice.WorkflowJsonSupport._
import cromwell.webservice.{EngineStatsActor, FailureResponse}
import mouse.boolean._
import org.scalatest.{AsyncFlatSpec, Matchers}
import spray.json._

import scala.concurrent.duration._

class CromwellApiServiceSpec extends AsyncFlatSpec with ScalatestRouteTest with Matchers {
  import CromwellApiServiceSpec._

  val akkaHttpService = new MockApiService()
  val version = "v1"

  implicit def default = RouteTestTimeout(5.seconds)

  "REST ENGINE /stats endpoint" should "no longer return 200 for stats" in {
    Get(s"/engine/$version/stats") ~>
      akkaHttpService.engineRoutes ~>
      check {
        status should be(StatusCodes.Forbidden)
        contentType should be(ContentTypes.`application/json`)
        responseAs[FailureResponse] should be(
          FailureResponse("fail", "The /stats endpoint is currently disabled.", None)
        )
      }
  }

  "REST ENGINE /version endpoint" should "return 200 for version" in {
    Get(s"/engine/$version/version") ~>
      akkaHttpService.engineRoutes ~>
      check {
        status should be(StatusCodes.OK)
        contentType should be(ContentTypes.`application/json`)
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
          contentType should be(ContentTypes.`application/json`)
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
          assertResult {
            s"""|{
                |  "status": "error",
                |  "message": "Couldn't abort $workflowId because no workflow with that ID is in progress"
                |}
                |""".stripMargin.trim
          } {
            responseAs[String]
          }
          assertResult(ContentTypes.`application/json`)(contentType)
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
          assertResult(ContentTypes.`application/json`)(contentType)
        }
    }

    it should "return 200 Aborted for abort of a workflow which a workflow is in OnHold state" in {
      Post(s"/workflows/$version/${CromwellApiServiceSpec.OnHoldWorkflowId}/abort") ~>
        akkaHttpService.workflowRoutes ~>
      check {
        assertResult(
          s"""{"id":"${CromwellApiServiceSpec.OnHoldWorkflowId.toString}","status":"Aborted"}""") {
          responseAs[String]
        }
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult(ContentTypes.`application/json`)(contentType)
      }
    }

  it should "return 200 Aborted for abort of a workflow which a workflow is in Submitted state" in {
    Post(s"/workflows/$version/${CromwellApiServiceSpec.SubmittedWorkflowId}/abort") ~>
      akkaHttpService.workflowRoutes ~>
      check {
        assertResult(
          s"""{"id":"${CromwellApiServiceSpec.SubmittedWorkflowId.toString}","status":"Aborted"}""") {
          responseAs[String]
        }
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult(ContentTypes.`application/json`)(contentType)
      }
  }

    it should "return 200 Aborting for abort of a known workflow id which is currently running" in {
      Post(s"/workflows/$version/${CromwellApiServiceSpec.AbortingWorkflowId}/abort") ~>
        akkaHttpService.workflowRoutes ~>
        check {
          assertResult(
            s"""{"id":"${CromwellApiServiceSpec.AbortingWorkflowId.toString}","status":"Aborting"}""") {
            responseAs[String]
          }
          assertResult(StatusCodes.OK) {
            status
          }
          assertResult(ContentTypes.`application/json`)(contentType)
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
             |  "message": "Error(s): Error while validating workflow url: unknown protocol: htpps",
             |  "status": "fail"
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
             |  "message": "Unrecognized workflow ID: ${CromwellApiServiceSpec.UnrecognizedWorkflowId.toString}",
             |  "status": "fail"
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
        contentType should be(ContentTypes.`application/json`)
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
          assertResult(ContentTypes.`application/json`)(contentType)
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
          assertResult {
            """|{
               |  "status": "fail",
               |  "message": "Error(s): Invalid workflow options provided: Unsupported key/value pair in WorkflowOptions: defaultRuntimeOptions -> {\"cpu\":1}"
               |}
               |""".stripMargin.trim
          } {
            responseAs[String]
          }
          assertResult(ContentTypes.`application/json`)(contentType)
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
          assertResult {
            """|{
               |  "status": "fail",
               |  "message": "Error(s): Invalid workflow options provided: Unexpected end-of-input at input index 28 (line 3, position 1), expected '}':\n\n^\n"
               |}
               |""".stripMargin.trim
          } {
            responseAs[String]
          }
          assertResult(ContentTypes.`application/json`)(contentType)
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
          assertResult(ContentTypes.`application/json`)(contentType)
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
          assertResult(ContentTypes.`application/json`)(contentType)
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
          assertResult(ContentTypes.`application/json`)(contentType)
        }
    }

    behavior of "REST API /timing endpoint"
    it should "return 200 with an HTML document for the timings route" in {
      Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/timing") ~>
        akkaHttpService.workflowRoutes ~>
        check {
          assertResult(StatusCodes.OK) { status }
          assertResult(`text/html(UTF-8)`) { contentType }
          assertResult("<html>") {
            responseAs[String].substring(0, 6)
          }
        }
    }

    it should "return 404 when unrecognized workflow id is submitted" in {
      Get(s"/workflows/$version/${CromwellApiServiceSpec.UnrecognizedWorkflowId}/timing") ~>
        akkaHttpService.workflowRoutes ~>
        check {
          assertResult(StatusCodes.NotFound) { status }
          assertResult(
            s"""{
               |  "message": "Unrecognized workflow ID: ${CromwellApiServiceSpec.UnrecognizedWorkflowId.toString}",
               |  "status": "fail"
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
               |  "message": "Invalid workflow ID: 'foo'.",
               |  "status": "fail"
               |}""".stripMargin) {
            responseAs[String].parseJson.prettyPrint
          }
          assertResult(ContentTypes.`application/json`)(contentType)
        }
    }
}

object CromwellApiServiceSpec {
  val ExistingWorkflowId = WorkflowId.fromString("c4c6339c-8cc9-47fb-acc5-b5cb8d2809f5")
  val AbortedWorkflowId = WorkflowId.fromString("0574111c-c7d3-4145-8190-7a7ed8e8324a")
  val UnrecognizedWorkflowId = WorkflowId.fromString("2bdd06cc-e794-46c8-a897-4c86cedb6a06")
  val OnHoldWorkflowId = WorkflowId.fromString("fe6dbaf6-e15c-4438-813c-479b35867142")
  val SubmittedWorkflowId = WorkflowId.fromString("1d4ec2e2-7c1c-407d-9be8-75b39ab72dfd")
  val RunningWorkflowId = WorkflowId.fromString("dcfea0ab-9a3e-4bc0-ab82-794a37c4d484")
  val AbortingWorkflowId = WorkflowId.fromString("2e3503f5-24f5-4a01-a4d1-bb1088bb5c1e")
  val SucceededWorkflowId = WorkflowId.fromString("0cb43b8c-0259-4a19-b7fe-921ced326738")
  val FailedWorkflowId = WorkflowId.fromString("df501790-cef5-4df7-9b48-8760533e3136")
  val SummarizedWorkflowId = WorkflowId.fromString("f0000000-0000-0000-0000-000000000000")
  val RecognizedWorkflowIds = Set(ExistingWorkflowId, AbortedWorkflowId, OnHoldWorkflowId, RunningWorkflowId, AbortingWorkflowId, SucceededWorkflowId, FailedWorkflowId, SummarizedWorkflowId)
  val SummarizedWorkflowIds = Set(SummarizedWorkflowId)

  class MockApiService()(implicit val system: ActorSystem) extends CromwellApiService {
    override def actorRefFactory = system

    override val materializer = ActorMaterializer()
    override val ec = system.dispatcher
    override val workflowStoreActor = actorRefFactory.actorOf(Props(new MockWorkflowStoreActor()))
    override val serviceRegistryActor = actorRefFactory.actorOf(Props(new MockServiceRegistryActor()))
    override val workflowManagerActor = actorRefFactory.actorOf(Props(new MockWorkflowManagerActor()))
  }

  object MockServiceRegistryActor {

    private def fullMetadataResponse(workflowId: WorkflowId) = {
      List(
        MetadataEvent(MetadataKey(workflowId, None, "testKey1a"), MetadataValue("myValue1a", MetadataString)),
        MetadataEvent(MetadataKey(workflowId, None, "testKey1b"), MetadataValue("myValue1b", MetadataString)),
        MetadataEvent(MetadataKey(workflowId, None, "testKey2a"), MetadataValue("myValue2a", MetadataString))
      )
    }

    def responseMetadataValues(workflowId: WorkflowId, withKeys: List[String], withoutKeys: List[String]): JsObject = {
      def keyFilter(keys: List[String])(m: MetadataEvent) = keys.exists(k => m.key.key.startsWith(k))
      val events = fullMetadataResponse(workflowId)
        .filter(m => withKeys.isEmpty || keyFilter(withKeys)(m))
        .filter(m => withoutKeys.isEmpty || !keyFilter(withoutKeys)(m))

      MetadataBuilderActor.workflowMetadataResponse(workflowId, events, includeCallsIfEmpty = false, Map.empty)
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

  class MockServiceRegistryActor extends Actor with StrictLogging {
    import MockServiceRegistryActor._

    override def receive = {
      case QueryForWorkflowsMatchingParameters(parameters) =>
        val labels: Option[Map[String, String]] = {
          parameters.contains(("additionalQueryResultFields", "labels")).option(
            Map("key1" -> "label1", "key2" -> "label2"))
        }

        val response = WorkflowQuerySuccess(WorkflowQueryResponse(List(WorkflowQueryResult(ExistingWorkflowId.toString,
          None, Some(WorkflowSucceeded.toString), None, None, None, labels, Option("pid"), Option("rid"))), 1), None)
        sender ! response
      case ValidateWorkflowIdInMetadata(id) =>
        if (RecognizedWorkflowIds.contains(id)) sender ! MetadataService.RecognizedWorkflowId
        else sender ! MetadataService.UnrecognizedWorkflowId
      case ValidateWorkflowIdInMetadataSummaries(id) =>
        if (SummarizedWorkflowIds.contains(id)) sender ! MetadataService.RecognizedWorkflowId
        else sender ! MetadataService.UnrecognizedWorkflowId
      case GetCurrentStatus =>
        sender ! StatusCheckResponse(
          ok = true,
          systems = Map(
            "Engine Database" -> SubsystemStatus(ok = true, messages = None)))
      case request @ GetStatus(id) =>
        val status = id match {
          case OnHoldWorkflowId => WorkflowOnHold
          case RunningWorkflowId => WorkflowRunning
          case AbortingWorkflowId => WorkflowAborting
          case AbortedWorkflowId => WorkflowAborted
          case SucceededWorkflowId => WorkflowSucceeded
          case FailedWorkflowId => WorkflowFailed
          case _ => WorkflowSubmitted
        }
        sender ! BuiltMetadataResponse(request, MetadataBuilderActor.processStatusResponse(id, status))
      case request @ GetLabels(id) =>
        sender ! BuiltMetadataResponse(request, MetadataBuilderActor.processLabelsResponse(id, Map("key1" -> "label1", "key2" -> "label2")))
      case request @ WorkflowOutputs(id) =>
        val event = Vector(MetadataEvent(MetadataKey(id, None, "outputs:test.hello.salutation"), MetadataValue("Hello foo!", MetadataString)))
        sender ! BuiltMetadataResponse(request, MetadataBuilderActor.processOutputsResponse(id, event))
      case request @ GetLogs(id) =>
        sender ! BuiltMetadataResponse(request, MetadataBuilderActor.workflowMetadataResponse(id, logsEvents(id), includeCallsIfEmpty = false, Map.empty))
      case request @ GetMetadataAction(MetadataQuery(id, _, _, withKeys, withoutKeys, _)) =>
        val withKeysList = withKeys.map(_.toList).getOrElse(List.empty)
        val withoutKeysList = withoutKeys.map(_.toList).getOrElse(List.empty)
        sender ! BuiltMetadataResponse(request, responseMetadataValues(id, withKeysList, withoutKeysList))
      case PutMetadataActionAndRespond(events, _, _) =>
        events.head.key.workflowId match {
          case CromwellApiServiceSpec.ExistingWorkflowId => sender ! MetadataWriteSuccess(events)
          case CromwellApiServiceSpec.SummarizedWorkflowId => sender ! MetadataWriteSuccess(events)
          case CromwellApiServiceSpec.AbortedWorkflowId => sender ! MetadataWriteFailure(new Exception("mock exception of db failure"), events)
          case WorkflowId(_) => throw new Exception("Something untoward happened, this situation is not believed to be possible at this time")
        }
      case DescribeRequest(sourceFiles) =>
        sourceFiles.workflowSource match {
          case Some("fail to describe") =>
            sender ! DescribeFailure("as requested, failing to describe")
          case Some("actor asplode") =>
            throw new Exception("asploding now!")
          case _ =>
            val readBack = List(
              "this is fake data from the mock SR actor",
              s"[reading back DescribeRequest contents] workflow hashcode: ${sourceFiles.workflowSource.map(_.hashCode)}",
              s"[reading back DescribeRequest contents] workflow url: ${sourceFiles.workflowUrl}",
              s"[reading back DescribeRequest contents] inputs: ${sourceFiles.inputsJson}",
              s"[reading back DescribeRequest contents] type: ${sourceFiles.workflowType}",
              s"[reading back DescribeRequest contents] version: ${sourceFiles.workflowTypeVersion}"
            )

            sender ! DescribeSuccess(description = WorkflowDescription(valid = true, errors = readBack, validWorkflow = true))
        }
      case _: InstrumentationServiceMessage => // Do nothing.
      case m => logger.error("Unexpected message received by MockServiceRegistryActor: {}", m)
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
          case AbortingWorkflowId => WorkflowAbortRequestedResponse(id)
          case OnHoldWorkflowId | SubmittedWorkflowId => WorkflowAbortedResponse(id)
          case UnrecognizedWorkflowId => WorkflowAbortFailureResponse(id, new WorkflowNotFoundException(s"Couldn't abort $id because no workflow with that ID is in progress"))
          case WorkflowId(_) => throw new Exception("Something untoward happened")
        }
        sender ! message
      case GetWorkflowStoreStats => sender ! Map(WorkflowRunning -> 5, WorkflowSubmitted -> 3, WorkflowAborting -> 2)
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
