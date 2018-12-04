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

import scala.concurrent.duration._


// N.B. this suite only tests the routing and initial validation, it uses the MockServiceRegistryActor
// to return fake results instead of going to WomtoolServiceActor
class WomtoolRouteSupportSpec extends AsyncFlatSpec with ScalatestRouteTest with Matchers {

  val akkaHttpService = new MockWomtoolRouteSupport()
  val version = "v1"

  implicit def default = RouteTestTimeout(5.seconds)

  behavior of "/describe endpoint"

  object TestEntities {
    val empty = HttpEntity.empty(ContentTypes.`application/json`)
  }

  it should "return OK for a correctly structured request" in {

    val workflow = Multipart.FormData.BodyPart("workflowSource", HttpEntity(MediaTypes.`application/json`, "this is not a WDL!"))
    val form = Multipart.FormData(workflow).toEntity()

    Post(s"/womtool/$version/describe", form) ~>
      akkaHttpService.womtoolRoutes ~>
      check {
        status should be(StatusCodes.OK)

        assertResult {
          DescribeResponse(valid = false, List("this is fake data", "from the mock SR actor"))
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
