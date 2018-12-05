package cromwell.webservice.routes

import akka.actor.{ActorSystem, Props}
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model._
import akka.http.scaladsl.testkit.{RouteTestTimeout, ScalatestRouteTest}
import akka.stream.ActorMaterializer
import cromwell.services.womtool.WomtoolServiceMessages.DescribeResponse
import cromwell.services.womtool.WomtoolServiceMessages.JsonSupport.describeResponseFormat
import cromwell.webservice.routes.CromwellApiServiceSpec.MockServiceRegistryActor
import cromwell.webservice.routes.WomtoolRouteSupportSpec.MockWomtoolRouteSupport
import org.scalatest.{AsyncFlatSpec, Matchers}
import spray.json._

import scala.concurrent.duration._


// N.B. this suite only tests the routing and initial validation, it uses the MockServiceRegistryActor
// to return fake results instead of going to WomtoolServiceActor
class WomtoolRouteSupportSpec extends AsyncFlatSpec with ScalatestRouteTest with Matchers {

  val akkaHttpService = new MockWomtoolRouteSupport()
  val version = "v1"

  implicit def default = RouteTestTimeout(5.seconds)

  behavior of "/describe endpoint"

  object BodyParts {
    val workflowSource  = Multipart.FormData.BodyPart("workflowSource", HttpEntity(MediaTypes.`application/json`, "This is not a WDL, but that's OK for this test of request routing."))
    val workflowUrl     = Multipart.FormData.BodyPart("workflowUrl", HttpEntity(MediaTypes.`application/json`,
      "https://raw.githubusercontent.com/broadinstitute/cromwell/develop/womtool/src/test/resources/validate/wdl_draft3/valid/callable_imports/my_workflow.wdl"))
    val workflowUrlNotFound = Multipart.FormData.BodyPart("workflowUrl", HttpEntity(MediaTypes.`application/json`,
      "https://raw.githubusercontent.com/broadinstitute/cromwell/develop/my_workflow"))
    val workflowUrlBadHost  = Multipart.FormData.BodyPart("workflowUrl", HttpEntity(MediaTypes.`application/json`, "https://zardoz.zardoz"))
    val workflowNotAUrl = Multipart.FormData.BodyPart("workflowUrl", HttpEntity(MediaTypes.`application/json`, "Zardoz"))
    val workflowInputs = Multipart.FormData.BodyPart("workflowInputs", HttpEntity(MediaTypes.`application/json`, "{\"a\":\"is for apple\"}"))
    val workflowType = Multipart.FormData.BodyPart("workflowType", HttpEntity(MediaTypes.`application/json`, "WDL"))
    val workflowVersion = Multipart.FormData.BodyPart("workflowTypeVersion", HttpEntity(MediaTypes.`application/json`, "1.0"))
  }

  it should "return OK for a workflow source request" in {
    Post(s"/womtool/$version/describe", Multipart.FormData(BodyParts.workflowSource).toEntity()) ~>
      akkaHttpService.womtoolRoutes ~>
      check {
        status should be(StatusCodes.OK)

        assertResult {
          DescribeResponse(valid = true, List("this is fake data", "from the mock SR actor", "workflow hashcode: 580529622, inputs: , type: None, version: None"))
        } { responseAs[DescribeResponse] }
      }
  }

  it should "return OK for a workflow URL request with a valid URL" in {
    Post(s"/womtool/$version/describe", Multipart.FormData(BodyParts.workflowUrl).toEntity()) ~>
      akkaHttpService.womtoolRoutes ~>
      check {
        status should be(StatusCodes.OK)

        assertResult {
          DescribeResponse(valid = true, List("this is fake data", "from the mock SR actor", "workflow hashcode: -1015985683, inputs: , type: None, version: None"))
        } { responseAs[DescribeResponse] }
      }
  }

  it should "return Bad Request when both URL and source provided" in {
    Post(s"/womtool/$version/describe", Multipart.FormData(BodyParts.workflowUrl, BodyParts.workflowSource).toEntity()) ~>
      akkaHttpService.womtoolRoutes ~>
      check {
        status should be(StatusCodes.BadRequest)

        assertResult(
          s"""{
             |  "status": "fail",
             |  "message": "Both workflow source and url can't be supplied"
             |}""".stripMargin) {
          responseAs[String]
        }
      }
  }

  it should "return Bad Request when neither URL or source provided" in {
    Post(s"/womtool/$version/describe", Multipart.FormData().toEntity()) ~>
      akkaHttpService.womtoolRoutes ~>
      check {
        status should be(StatusCodes.BadRequest)

        assertResult(
          s"""{
             |  "status": "fail",
             |  "message": "Either workflow source or url has to be supplied"
             |}""".stripMargin) {
          responseAs[String]
        }
      }
  }

  it should "return Bad Request when the workflow URL is a 404" in {
    Post(s"/womtool/$version/describe", Multipart.FormData(BodyParts.workflowUrlNotFound).toEntity()) ~>
      akkaHttpService.womtoolRoutes ~>
      check {
        status should be(StatusCodes.BadRequest)

        assertResult(
          s"""{
             |  "status": "fail",
             |  "message": "Failed to resolve 'https://raw.githubusercontent.com/broadinstitute/cromwell/develop/my_workflow' using resolver: 'http importer (no 'relative-to' origin)' (reason 1 of 1): Failed to download https://raw.githubusercontent.com/broadinstitute/cromwell/develop/my_workflow (reason 1 of 1): 404: Not Found\\n"
             |}""".stripMargin) {
          responseAs[String]
        }
      }
  }

  it should "return Bad Request when the workflow URL's host can't be resolved" in {
    Post(s"/womtool/$version/describe", Multipart.FormData(BodyParts.workflowUrlBadHost).toEntity()) ~>
      akkaHttpService.womtoolRoutes ~>
      check {
        status should be(StatusCodes.BadRequest)

        val responseJson: JsObject = responseAs[String].parseJson.asJsObject

        responseJson.fields("status") shouldBe JsString("fail")

        // The error is from the OS network stack and differs between Mac and Linux
        responseJson.fields("message") should (
          be(JsString("Failed to resolve 'https://zardoz.zardoz' using resolver: 'http importer (no 'relative-to' origin)' (reason 1 of 1): Failed to download https://zardoz.zardoz (reason 1 of 1): HTTP resolver with headers had an unexpected error (zardoz.zardoz: Name or service not known)"))
          or
          be(JsString("Failed to resolve 'https://zardoz.zardoz' using resolver: 'http importer (no 'relative-to' origin)' (reason 1 of 1): Failed to download https://zardoz.zardoz (reason 1 of 1): HTTP resolver with headers had an unexpected error (zardoz.zardoz: nodename nor servname provided, or not known)"))
        )
      }
  }

  it should "return Bad Request when the workflow URL is not a URL" in {
    Post(s"/womtool/$version/describe", Multipart.FormData(BodyParts.workflowNotAUrl).toEntity()) ~>
      akkaHttpService.womtoolRoutes ~>
      check {
        status should be(StatusCodes.BadRequest)

        assertResult(
          s"""{
             |  "status": "fail",
             |  "message": "Failed to resolve 'Zardoz' using resolver: 'http importer (no 'relative-to' origin)' (reason 1 of 1): Relative path"
             |}""".stripMargin) {
          responseAs[String]
        }
      }
  }

  it should "include inputs, workflow type, and workflow version in the WorkflowSourceFilesCollection" in {
    Post(s"/womtool/$version/describe", Multipart.FormData(BodyParts.workflowSource, BodyParts.workflowInputs, BodyParts.workflowType, BodyParts.workflowVersion).toEntity()) ~>
      akkaHttpService.womtoolRoutes ~>
      check {
        status should be(StatusCodes.OK)

        assertResult {
          DescribeResponse(valid = true, List("this is fake data", "from the mock SR actor", "workflow hashcode: 580529622, inputs: {\"a\":\"is for apple\"}, type: Some(WDL), version: Some(1.0)"))
        } { responseAs[DescribeResponse] }
      }
  }
}

object WomtoolRouteSupportSpec {
  class MockWomtoolRouteSupport()(implicit val system: ActorSystem, routeTestTimeout: RouteTestTimeout) extends WomtoolRouteSupport {
    override def actorRefFactory = system
    override val ec = system.dispatcher
    override val timeout = routeTestTimeout.duration
    override val serviceRegistryActor = actorRefFactory.actorOf(Props(new MockServiceRegistryActor()))
    override implicit val materializer = ActorMaterializer()
  }
}
