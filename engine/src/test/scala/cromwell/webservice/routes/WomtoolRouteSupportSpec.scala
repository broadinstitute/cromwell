package cromwell.webservice.routes

import akka.actor.{ActorSystem, Props}
import akka.http.scaladsl.model._
import akka.http.scaladsl.testkit.{RouteTestTimeout, ScalatestRouteTest}
import akka.http.scaladsl.unmarshalling.Unmarshaller._
import akka.stream.ActorMaterializer
import akka.testkit._
import cromwell.services.womtool.models.WorkflowDescription
import cromwell.services.womtool.models.WorkflowDescription.workflowDescriptionDecoder
import cromwell.webservice.routes.CromwellApiServiceSpec.MockServiceRegistryActor
import cromwell.webservice.routes.WomtoolRouteSupportSpec.MockWomtoolRouteSupport
import de.heikoseeberger.akkahttpcirce.ErrorAccumulatingCirceSupport._
import org.scalatest.{AsyncFlatSpec, Matchers}

import scala.concurrent.duration._


// N.B. this suite only tests the routing and creation of the WorkflowSourceFilesCollection, it uses the MockServiceRegistryActor
// to return fake results instead of going to a real WomtoolServiceActor
class WomtoolRouteSupportSpec extends AsyncFlatSpec with ScalatestRouteTest with Matchers {

  val akkaHttpService = new MockWomtoolRouteSupport()
  val version = "v1"

  implicit def default = RouteTestTimeout(10.seconds.dilated)

  behavior of "/describe endpoint"

  object BodyParts {
    val workflowSource  = Multipart.FormData.BodyPart("workflowSource", HttpEntity(ContentTypes.`text/plain(UTF-8)`, "This is not a WDL, but that's OK for this test of request routing."))
    val workflowUrl     = Multipart.FormData.BodyPart("workflowUrl", HttpEntity(MediaTypes.`application/json`,
      "https://raw.githubusercontent.com/broadinstitute/cromwell/develop/womtool/src/test/resources/validate/wdl_draft3/valid/callable_imports/my_workflow.wdl"))
    val workflowInputs = Multipart.FormData.BodyPart("workflowInputs", HttpEntity(MediaTypes.`application/json`, "{\"a\":\"is for apple\"}"))
    val workflowType = Multipart.FormData.BodyPart("workflowType", HttpEntity(ContentTypes.`text/plain(UTF-8)`, "WDL"))
    val workflowVersion = Multipart.FormData.BodyPart("workflowTypeVersion", HttpEntity(ContentTypes.`text/plain(UTF-8)`, "1.0"))

    val workflowSourceTriggerDescribeFailure =
      Multipart.FormData.BodyPart("workflowSource", HttpEntity(ContentTypes.`text/plain(UTF-8)`, "fail to describe"))
    val workflowSourceTriggerActorException =
      Multipart.FormData.BodyPart("workflowSource", HttpEntity(ContentTypes.`text/plain(UTF-8)`, "actor asplode"))
  }

  it should "return Bad Request if the actor returns DescribeFailure" in {
    Post(s"/womtool/$version/describe", Multipart.FormData(BodyParts.workflowSourceTriggerDescribeFailure).toEntity()) ~>
      akkaHttpService.womtoolRoutes ~>
      check {
        status should be(StatusCodes.BadRequest)

        assertResult(
          s"""{
             |  "status": "fail",
             |  "message": "as requested, failing to describe"
             |}""".stripMargin) {
          responseAs[String]
        }
      }
  }

  it should "return Internal Server Error if the actor explodes" in {
    Post(s"/womtool/$version/describe", Multipart.FormData(BodyParts.workflowSourceTriggerActorException).toEntity()) ~>
      akkaHttpService.womtoolRoutes ~>
      check {
        status should be(StatusCodes.InternalServerError)

        responseAs[String] should include("Ask timed out on")
      }
  }

  it should "return OK for a workflow source request" in {
    Post(s"/womtool/$version/describe", Multipart.FormData(BodyParts.workflowSource).toEntity()) ~>
      akkaHttpService.womtoolRoutes ~>
      check {
        status should be(StatusCodes.OK)

        assertResult {
          WorkflowDescription(valid = true,
            errors = List(
              "this is fake data from the mock SR actor",
              "[reading back DescribeRequest contents] workflow hashcode: Some(580529622)",
              "[reading back DescribeRequest contents] workflow url: None",
              "[reading back DescribeRequest contents] inputs: ",
              "[reading back DescribeRequest contents] type: None",
              "[reading back DescribeRequest contents] version: None"
            ),
            validWorkflow = true
          )
        } { responseAs[WorkflowDescription] }
      }
  }

  it should "return OK for a workflow URL request with a valid URL" in {
    Post(s"/womtool/$version/describe", Multipart.FormData(BodyParts.workflowUrl).toEntity()) ~>
      akkaHttpService.womtoolRoutes ~>
      check {
        status should be(StatusCodes.OK)

        assertResult {
          WorkflowDescription(valid = true,
            errors = List(
              "this is fake data from the mock SR actor",
              "[reading back DescribeRequest contents] workflow hashcode: None",
              "[reading back DescribeRequest contents] workflow url: Some(https://raw.githubusercontent.com/broadinstitute/cromwell/develop/womtool/src/test/resources/validate/wdl_draft3/valid/callable_imports/my_workflow.wdl)",
              "[reading back DescribeRequest contents] inputs: ",
              "[reading back DescribeRequest contents] type: None",
              "[reading back DescribeRequest contents] version: None"
            ),
            validWorkflow = true
          )
        } { responseAs[WorkflowDescription] }
      }
  }

  it should "include inputs, workflow type, and workflow version in the WorkflowSourceFilesCollection" in {
    Post(s"/womtool/$version/describe", Multipart.FormData(BodyParts.workflowSource, BodyParts.workflowInputs, BodyParts.workflowType, BodyParts.workflowVersion).toEntity()) ~>
      akkaHttpService.womtoolRoutes ~>
      check {
        status should be(StatusCodes.OK)

        assertResult {
          WorkflowDescription(valid = true,
            errors = List(
              "this is fake data from the mock SR actor",
              "[reading back DescribeRequest contents] workflow hashcode: Some(580529622)",
              "[reading back DescribeRequest contents] workflow url: None",
              "[reading back DescribeRequest contents] inputs: {\"a\":\"is for apple\"}",
              "[reading back DescribeRequest contents] type: Some(WDL)",
              "[reading back DescribeRequest contents] version: Some(1.0)"
            ),
            validWorkflow = true
          )
        } { responseAs[WorkflowDescription] }
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
