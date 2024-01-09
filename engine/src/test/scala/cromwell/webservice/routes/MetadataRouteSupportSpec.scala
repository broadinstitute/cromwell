package cromwell.webservice.routes

import akka.actor.{ActorSystem, Props}
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model.headers.{`Accept-Encoding`, HttpEncodings}
import akka.http.scaladsl.model.{ContentTypes, HttpEntity, StatusCodes}
import akka.http.scaladsl.server.Route
import akka.http.scaladsl.testkit.{RouteTestTimeout, ScalatestRouteTest}
import cromwell.core.WorkflowMetadataKeys
import cromwell.webservice.routes.CromwellApiServiceSpec.MockServiceRegistryActor
import cromwell.webservice.routes.MetadataRouteSupportSpec.MockMetadataRouteSupport
import org.scalatest.flatspec.AsyncFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json.DefaultJsonProtocol._
import spray.json._

import scala.concurrent.duration._

class MetadataRouteSupportSpec extends AsyncFlatSpec with ScalatestRouteTest with Matchers {
  val akkaHttpService = new MockMetadataRouteSupport()

  val version = "v1"

  implicit def routeTestTimeout = RouteTestTimeout(5.seconds)

  behavior of "REST API /status endpoint"
  it should "return 200 for get of a known workflow id" in {
    val workflowId = CromwellApiServiceSpec.ExistingWorkflowId

    Get(s"/workflows/$version/$workflowId/status") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status should be(StatusCodes.OK)
        // Along w/ checking value, ensure it is valid JSON despite the requested content type
        responseAs[JsObject].fields(WorkflowMetadataKeys.Status) should be(JsString("Submitted"))
        contentType should be(ContentTypes.`application/json`)
      }
  }

  it should "return 404 for get of unknown workflow" in {
    val workflowId = CromwellApiServiceSpec.UnrecognizedWorkflowId
    Get(s"/workflows/$version/$workflowId/status") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
        assertResult(
          s"""|{
              |  "status": "fail",
              |  "message": "Unrecognized workflow ID: ${CromwellApiServiceSpec.UnrecognizedWorkflowId}"
              |}
              |""".stripMargin.trim
        ) {
          responseAs[String]
        }
        assertResult(ContentTypes.`application/json`)(contentType)
      }
  }

  it should "return 400 for get of a malformed workflow id's status" in {
    Get(s"/workflows/$version/foobar/status") ~>
      akkaHttpService.metadataRoutes ~>
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

  behavior of "REST API /outputs endpoint"
  it should "return 200 with GET of outputs on successful execution of workflow" in {
    Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/outputs") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status should be(StatusCodes.OK)
        responseAs[JsObject].fields.keys should contain allOf (WorkflowMetadataKeys.Id, WorkflowMetadataKeys.Outputs)
        contentType should be(ContentTypes.`application/json`)
      }
  }

  it should "return 404 with outputs on unknown workflow" in {
    Get(s"/workflows/$version/${CromwellApiServiceSpec.UnrecognizedWorkflowId}/outputs") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
        assertResult(ContentTypes.`application/json`)(contentType)
      }
  }

  it should "return 405 with POST of outputs on successful execution of workflow" in {
    Post(s"/workflows/$version/${CromwellApiServiceSpec.UnrecognizedWorkflowId}/outputs") ~>
      Route.seal(akkaHttpService.metadataRoutes) ~>
      check {
        assertResult(StatusCodes.MethodNotAllowed) {
          status
        }
        assertResult("HTTP method not allowed, supported methods: GET")(responseAs[String])
        assertResult(ContentTypes.`text/plain(UTF-8)`)(contentType)
      }
  }

  it should "return 200 with GET of outputs of archived workflow (not deleted)" in {
    Get(s"/workflows/$version/${CromwellApiServiceSpec.ArchivedWorkflowId}/outputs") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status should be(StatusCodes.OK)
        responseAs[JsObject].fields.keys should contain allOf (WorkflowMetadataKeys.Id, WorkflowMetadataKeys.Outputs)
        contentType should be(ContentTypes.`application/json`)
      }
  }

  def validateArchivedMetadataResponseMessage(responseJson: JsObject,
                                              includeAvailabilityMessage: Boolean,
                                              includeLabelsMessage: Boolean
  ) = {
    val responseMessage = responseJson.fields("message").asInstanceOf[JsString].value
    val expectedSuffix =
      (if (includeAvailabilityMessage)
         " It is available in the archive bucket, or via a support request in the case of a managed instance."
       else "") +
        (if (includeLabelsMessage)
           " As a result, new labels can't be added or existing labels can't be updated for this workflow."
         else "")

    responseMessage should startWith(
      "Cromwell has archived this workflow's metadata " +
        "according to the lifecycle policy. The workflow completed at "
    )
    // The missing middle of the message looks like "The workflow completed at x timestamp, which was y milliseconds ago."
    responseMessage should endWith(s"ago.${expectedSuffix}")
  }

  it should "return 200 with ArchivedAndDeleted status for workflow whose metadata is deleted" in {
    Get(s"/workflows/$version/${CromwellApiServiceSpec.ArchivedAndDeletedWorkflowId}/outputs") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status should be(StatusCodes.OK)

        val responseJson = responseAs[JsObject]
        responseJson.fields.keys should contain allOf (WorkflowMetadataKeys.Id, WorkflowMetadataKeys.MetadataArchiveStatus, WorkflowMetadataKeys.Message)
        responseJson.fields("metadataArchiveStatus").asInstanceOf[JsString].value shouldBe "ArchivedAndDeleted"
        validateArchivedMetadataResponseMessage(responseJson,
                                                includeAvailabilityMessage = true,
                                                includeLabelsMessage = false
        )
        contentType should be(ContentTypes.`application/json`)
      }
  }

  behavior of "REST API /logs endpoint"
  it should "return 200 with paths to stdout/stderr/backend log" in {
    Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/logs") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status should be(StatusCodes.OK)

        val call =
          responseAs[JsObject].fields("calls").convertTo[JsObject].fields("mycall").convertTo[Seq[JsObject]].head
        call.fields("stdout") should be(JsString("stdout.txt"))
        call.fields("stderr") should be(JsString("stderr.txt"))
        call.fields("stdout") should be(JsString("stdout.txt"))
        call.fields("backendLogs").convertTo[JsObject].fields("log") should be(JsString("backend.log"))
      }
  }

  it should "return 404 with logs on unknown workflow" in {
    Get(s"/workflows/$version/${CromwellApiServiceSpec.UnrecognizedWorkflowId}/logs") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
      }
  }

  it should "return 200 with paths to stdout/stderr/backend log for archived workflow (not deleted)" in {
    Get(s"/workflows/$version/${CromwellApiServiceSpec.ArchivedWorkflowId}/logs") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status should be(StatusCodes.OK)

        val call =
          responseAs[JsObject].fields("calls").convertTo[JsObject].fields("mycall").convertTo[Seq[JsObject]].head
        call.fields("stdout") should be(JsString("stdout.txt"))
        call.fields("stderr") should be(JsString("stderr.txt"))
        call.fields("stdout") should be(JsString("stdout.txt"))
        call.fields("backendLogs").convertTo[JsObject].fields("log") should be(JsString("backend.log"))
      }
  }

  it should "return 200 with ArchivedAndDeleted status for workflow whose metadata is deleted" in {
    Get(s"/workflows/$version/${CromwellApiServiceSpec.ArchivedAndDeletedWorkflowId}/logs") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status should be(StatusCodes.OK)

        val responseJson = responseAs[JsObject]
        responseJson.fields.keys should contain allOf (WorkflowMetadataKeys.Id, WorkflowMetadataKeys.MetadataArchiveStatus, WorkflowMetadataKeys.Message)
        responseJson.fields("metadataArchiveStatus").asInstanceOf[JsString].value shouldBe "ArchivedAndDeleted"
        validateArchivedMetadataResponseMessage(responseJson,
                                                includeAvailabilityMessage = true,
                                                includeLabelsMessage = false
        )
      }
  }

  behavior of "REST API /metadata endpoint"
  it should "return with full metadata from the metadata route" in {
    Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/metadata") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status should be(StatusCodes.OK)
        val result = responseAs[JsObject]
        result.fields.keys should contain allOf ("testKey1a", "testKey1b", "testKey2a")
        result.fields.keys shouldNot contain("testKey3")
        result.fields("testKey1a") should be(JsString("myValue1a"))
        result.fields("testKey1b") should be(JsString("myValue1b"))
        result.fields("testKey2a") should be(JsString("myValue2a"))
      }
  }

  it should "return with full metadata from the metadata route for workflow id which exists only in summary table" in {
    Get(s"/workflows/$version/${CromwellApiServiceSpec.WorkflowIdExistingOnlyInSummaryTable}/metadata") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status should be(StatusCodes.OK)
        val result = responseAs[JsObject]
        result.fields.keys should contain allOf ("testKey1a", "testKey1b", "testKey2a")
        result.fields.keys shouldNot contain("testKey3")
        result.fields("testKey1a") should be(JsString("myValue1a"))
        result.fields("testKey1b") should be(JsString("myValue1b"))
        result.fields("testKey2a") should be(JsString("myValue2a"))
      }
  }

  it should "return with gzip encoding when requested" in {
    Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/metadata")
      .addHeader(`Accept-Encoding`(HttpEncodings.gzip)) ~>
      akkaHttpService.metadataRoutes ~>
      check {
        response.headers.find(_.name == "Content-Encoding").get.value should be("gzip")
      }
  }

  it should "not return with gzip encoding when not requested" in {
    Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/metadata") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        response.headers.find(_.name == "Content-Encoding") shouldBe None
      }
  }

  it should "return with included metadata from the metadata route" in {
    Get(
      s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/metadata?includeKey=testKey1&includeKey=testKey2a"
    ) ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status should be(StatusCodes.OK)
        val result = responseAs[JsObject]
        result.fields.keys should contain allOf ("testKey1a", "testKey1b", "testKey2a")
        result.fields.keys should contain noneOf ("testKey2b", "testKey3")
        result.fields("testKey1a") should be(JsString("myValue1a"))
        result.fields("testKey1b") should be(JsString("myValue1b"))
        result.fields("testKey2a") should be(JsString("myValue2a"))
      }
  }

  it should "return with excluded metadata from the metadata route" in {
    Get(
      s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/metadata?excludeKey=testKey2&excludeKey=testKey3"
    ) ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status should be(StatusCodes.OK)
        val result = responseAs[JsObject]
        result.fields.keys should contain allOf ("testKey1a", "testKey1b")
        result.fields.keys should contain noneOf ("testKey2a", "testKey3")
        result.fields("testKey1a") should be(JsString("myValue1a"))
        result.fields("testKey1b") should be(JsString("myValue1b"))
      }
  }

  it should "correctly include and exclude metadata keys in workflow details requests" in {
    Get(
      s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/metadata?includeKey=testKey1&excludeKey=testKey1a"
    ) ~>
      akkaHttpService.metadataRoutes ~>
      check {
        val r = responseAs[String]
        withClue(s"From response $r") {
          status should be(StatusCodes.OK)
          val result = responseAs[JsObject]
          result.fields.keys should contain allElementsOf (List("testKey1b"))
          result.fields.keys should contain noneOf ("testKey1a", "testKey2")
          result.fields("testKey1b") should be(JsString("myValue1b"))
        }
      }
  }

  it should "return with full metadata from the metadata route for archived workflow (not deleted)" in {
    Get(s"/workflows/$version/${CromwellApiServiceSpec.ArchivedWorkflowId}/metadata") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status should be(StatusCodes.OK)
        val result = responseAs[JsObject]
        result.fields.keys should contain allOf ("testKey1a", "testKey1b", "testKey2a")
        result.fields.keys shouldNot contain("testKey3")
        result.fields("testKey1a") should be(JsString("myValue1a"))
        result.fields("testKey1b") should be(JsString("myValue1b"))
        result.fields("testKey2a") should be(JsString("myValue2a"))
      }
  }

  it should "return 200 with ArchivedAndDeleted status for workflow whose metadata is deleted" in {
    Get(s"/workflows/$version/${CromwellApiServiceSpec.ArchivedAndDeletedWorkflowId}/metadata") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status should be(StatusCodes.OK)

        val responseJson = responseAs[JsObject]
        responseJson.fields.keys should contain allOf (WorkflowMetadataKeys.Id, WorkflowMetadataKeys.MetadataArchiveStatus, WorkflowMetadataKeys.Message)
        responseJson.fields("metadataArchiveStatus").asInstanceOf[JsString].value shouldBe "ArchivedAndDeleted"
        validateArchivedMetadataResponseMessage(responseJson,
                                                includeAvailabilityMessage = true,
                                                includeLabelsMessage = false
        )
      }
  }

  it should "return 200 with metadata for failed tasks on a workflow" in {
    Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/metadata/failed-jobs") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status should be(StatusCodes.OK)
        val result = responseAs[JsObject]
        result.fields.keys should contain allOf ("testKey1a", "testKey1b", "testKey2a")
        result.fields.keys shouldNot contain("testKey3")
        result.fields("testKey1a") should be(JsString("myValue1a"))
        result.fields("testKey1b") should be(JsString("myValue1b"))
        result.fields("testKey2a") should be(JsString("myValue2a"))
      }
  }

  behavior of "REST API /query GET endpoint"
  it should "return good results for a good query" in {
    Get(s"/workflows/$version/query?status=Succeeded&id=${CromwellApiServiceSpec.ExistingWorkflowId}") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status should be(StatusCodes.OK)
        contentType should be(ContentTypes.`application/json`)
        val results = responseAs[JsObject].fields("results").convertTo[Seq[JsObject]]
        results.head.fields("id") should be(JsString(CromwellApiServiceSpec.ExistingWorkflowId.toString))
        results.head.fields("status") should be(JsString("Succeeded"))
      }
  }

  it should "return labels if specified in additionalQueryResultFields param" in {
    Get(
      s"/workflows/$version/query?additionalQueryResultFields=labels&id=${CromwellApiServiceSpec.ExistingWorkflowId}"
    ) ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status should be(StatusCodes.OK)
        contentType should be(ContentTypes.`application/json`)
        val results = responseAs[JsObject].fields("results").convertTo[Seq[JsObject]]
        val fields = results.head.fields
        fields("id") should be(JsString(CromwellApiServiceSpec.ExistingWorkflowId.toString))
        fields(WorkflowMetadataKeys.Labels).asJsObject.fields("key1") should be(JsString("label1"))
        fields(WorkflowMetadataKeys.Labels).asJsObject.fields("key2") should be(JsString("label2"))
      }
  }

  it should "return parentWorkflowId if specified in additionalQueryResultFields param" in {
    Get(
      s"/workflows/$version/query?additionalQueryResultFields=parentWorkflowId&id=${CromwellApiServiceSpec.ExistingWorkflowId}"
    ) ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status should be(StatusCodes.OK)
        contentType should be(ContentTypes.`application/json`)
        val results = responseAs[JsObject].fields("results").convertTo[Seq[JsObject]]
        val fields = results.head.fields
        fields("id") should be(JsString(CromwellApiServiceSpec.ExistingWorkflowId.toString))
        fields(WorkflowMetadataKeys.ParentWorkflowId) should be(JsString("pid"))
      }
  }

  behavior of "REST API /query POST endpoint"
  it should "return good results for a good query map body" in {
    Post(s"/workflows/$version/query", HttpEntity(ContentTypes.`application/json`, """[{"status":"Succeeded"}]""")) ~>
      akkaHttpService.metadataRoutes ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult(true) {
          entityAs[String].contains("\"status\":\"Succeeded\"")
        }
      }
  }

  it should "return labels if specified in additionalQueryResultFields param" in {
    Post(s"/workflows/$version/query",
         HttpEntity(ContentTypes.`application/json`, """[{"additionalQueryResultFields":"labels"}]""")
    ) ~>
      akkaHttpService.metadataRoutes ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        val results = responseAs[JsObject].fields("results").convertTo[Seq[JsObject]]
        val fields = results.head.fields
        fields(WorkflowMetadataKeys.Labels).asJsObject.fields("key1") should be(JsString("label1"))
        fields(WorkflowMetadataKeys.Labels).asJsObject.fields("key2") should be(JsString("label2"))
      }
  }

  it should "return parentWorkflowId if specified in additionalQueryResultFields param" in {
    Post(s"/workflows/$version/query",
         HttpEntity(ContentTypes.`application/json`, """[{"additionalQueryResultFields":"parentWorkflowId"}]""")
    ) ~>
      akkaHttpService.metadataRoutes ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult(true) {
          entityAs[String].contains("\"parentWorkflowId\":\"pid\"")
        }
      }
  }

  it should "include totalResultCount in workflow query response" in {
    Post(s"/workflows/$version/query", HttpEntity(ContentTypes.`application/json`, """[{"status":"Succeeded"}]""")) ~>
      akkaHttpService.metadataRoutes ~>
      check {
        val result = responseAs[JsObject]
        status should be(StatusCodes.OK)
        result.fields("totalResultsCount") should be(JsNumber("1"))
      }
  }

  behavior of "REST API /labels GET endpoint"
  it should "return labels for a workflow ID" in {
    Get(s"/workflows/$version/${CromwellApiServiceSpec.ExistingWorkflowId}/labels") ~>
      akkaHttpService.metadataRoutes ~>
      check {
        val result = responseAs[JsObject]
        status should be(StatusCodes.OK)
        result.fields(WorkflowMetadataKeys.Labels).asJsObject.fields("key1") should be(JsString("label1"))
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

    val workflowId = CromwellApiServiceSpec.SummarizedWorkflowId

    Patch(s"/workflows/$version/$workflowId/labels", HttpEntity(ContentTypes.`application/json`, validLabelsJson)) ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status shouldBe StatusCodes.OK
        val actualResult = responseAs[JsObject]
        val expectedResults =
          s"""
             |{
             |  "id": "$workflowId",
             |  "labels": {
             |    "key1": "label1",
             |    "key2": "label2",
             |    "label-key-1":"label-value-1",
             |    "label-key-2":"label-value-2"
             |  }
             |}
            """.stripMargin.parseJson

        actualResult shouldBe expectedResults
      }
  }

  it should "not update labels for workflow whose metadata is archived (not deleted)" in {
    val validLabelsJson =
      """
        |{
        |  "label-key-1":"label-value-1",
        |  "label-key-2":"label-value-2"
        |}
      """.stripMargin

    val workflowId = CromwellApiServiceSpec.ArchivedWorkflowId

    Patch(s"/workflows/$version/$workflowId/labels", HttpEntity(ContentTypes.`application/json`, validLabelsJson)) ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status shouldBe StatusCodes.BadRequest
        val actualResult = responseAs[JsObject]
        actualResult.fields.keys should contain allOf (WorkflowMetadataKeys.Id, WorkflowMetadataKeys.MetadataArchiveStatus, WorkflowMetadataKeys.Message)
        actualResult.fields("metadataArchiveStatus").asInstanceOf[JsString].value shouldBe "Archived"
        validateArchivedMetadataResponseMessage(actualResult,
                                                includeAvailabilityMessage = false,
                                                includeLabelsMessage = true
        )
      }
  }

  it should "not update labels for workflow whose metadata is archived and deleted" in {
    val validLabelsJson =
      """
        |{
        |  "label-key-1":"label-value-1",
        |  "label-key-2":"label-value-2"
        |}
      """.stripMargin

    val workflowId = CromwellApiServiceSpec.ArchivedAndDeletedWorkflowId

    Patch(s"/workflows/$version/$workflowId/labels", HttpEntity(ContentTypes.`application/json`, validLabelsJson)) ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status shouldBe StatusCodes.BadRequest
        val actualResult = responseAs[JsObject]
        actualResult.fields.keys should contain allOf (WorkflowMetadataKeys.Id, WorkflowMetadataKeys.MetadataArchiveStatus, WorkflowMetadataKeys.Message)
        actualResult.fields("metadataArchiveStatus").asInstanceOf[JsString].value shouldBe "ArchivedAndDeleted"
        validateArchivedMetadataResponseMessage(actualResult,
                                                includeAvailabilityMessage = true,
                                                includeLabelsMessage = true
        )
      }
  }

  it should "fail when assigning valid labels to an unsummarized workflow ID" in {
    val validLabelsJson =
      """
        |{
        |  "label-key-1":"label-value-1",
        |  "label-key-2":"label-value-2"
        |}
      """.stripMargin

    val unsummarizedId = CromwellApiServiceSpec.ExistingWorkflowId
    Patch(s"/workflows/$version/$unsummarizedId/labels",
          HttpEntity(ContentTypes.`application/json`, validLabelsJson)
    ) ~>
      akkaHttpService.metadataRoutes ~>
      check {
        status shouldBe StatusCodes.NotFound
      }

  }
}

object MetadataRouteSupportSpec {
  class MockMetadataRouteSupport()(implicit val system: ActorSystem, routeTestTimeout: RouteTestTimeout)
      extends MetadataRouteSupport {
    override def actorRefFactory = system
    override val ec = system.dispatcher
    override val timeout = routeTestTimeout.duration
    override val serviceRegistryActor = actorRefFactory.actorOf(Props(new MockServiceRegistryActor()))
  }
}
