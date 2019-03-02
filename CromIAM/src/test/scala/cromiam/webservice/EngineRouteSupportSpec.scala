package cromiam.webservice

import akka.http.scaladsl.model.ContentTypes
import akka.http.scaladsl.model.StatusCodes._
import akka.http.scaladsl.model.headers.Authorization
import akka.http.scaladsl.testkit.ScalatestRouteTest
import cromiam.server.status.{MockStatusService, StatusService}
import org.scalatest.{FlatSpec, Matchers}

class EngineRouteSupportSpec extends FlatSpec with Matchers with ScalatestRouteTest with EngineRouteSupport  {
  override val cromwellClient = new MockCromwellClient()
  val samClient = new MockSamClient()
  override val statusService: StatusService = new MockStatusService(() => Map.empty)

  val version = "v1"

  val versionRoutePath: String = s"/engine/$version/version"
  val statsRoutePath: String = s"/engine/$version/stats"
  val engineStatusRoutePath: String = s"/engine/$version/status"

  val emptyAuthHeader: Authorization = new Authorization(null)

  "Stats endpoint" should "return Forbidden as it is disabled" in {
    Get(statsRoutePath).withHeaders(emptyAuthHeader) ~> engineRoutes ~> check {
      status shouldBe Forbidden
      response.entity shouldBe EngineRouteSupport.CromIamStatsForbidden.entity
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  "Version endpoint" should "return 200 without authorization" in {
    Get(versionRoutePath) ~> engineRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe "Response from Cromwell"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  "Status endpoint" should "return 200 and statuses of it's subsystems" in {
    Get(engineStatusRoutePath).withHeaders(emptyAuthHeader) ~> engineRoutes ~> check {
      status shouldBe OK
      assertResult("""{"ok":true,"systems":{"Cromwell":{"ok":true},"Sam":{"ok":true}}}""") {
        responseAs[String]
      }
      assertResult(ContentTypes.`application/json`)(contentType)
    }
  }
}
