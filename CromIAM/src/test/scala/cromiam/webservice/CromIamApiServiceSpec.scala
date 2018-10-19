package cromiam.webservice

import akka.event.NoLogging
import akka.http.scaladsl.model.StatusCodes._
import akka.http.scaladsl.model.headers.Authorization
import akka.http.scaladsl.testkit.ScalatestRouteTest
import cromiam.server.status.StatusService
import org.scalatest._

class CromIamApiServiceSpec extends FlatSpec with Matchers with CromIamApiService with ScalatestRouteTest {
  override def testConfigSource = "akka.loglevel = DEBUG"
  val log = NoLogging
  // We need an auth object, but it's never used, so...
  val emptyAuthObject = new Authorization(null)

  override def configuration = throw new NotImplementedError("This spec shouldn't need to access the real interface/port")

  override lazy val cromwellClient = new MockCromwellClient()

  override val statusService: StatusService = new StatusService(() => Map.empty)

  "Stats endpoint" should "be forbidden" in {
    Get("/engine/v1/stats").withHeaders(new Authorization(null)) ~> allRoutes ~> check {
      status shouldBe Forbidden
      response.entity shouldBe EngineRouteSupport.CromIamStatsForbidden.entity
    }
  }
}
