package cromiam.webservice

import akka.actor.ActorSystem
import akka.event.NoLogging
import akka.http.scaladsl.model.{HttpRequest, HttpResponse}
import akka.http.scaladsl.model.StatusCodes._
import akka.http.scaladsl.model.headers.Authorization
import akka.http.scaladsl.testkit.ScalatestRouteTest
import akka.stream.ActorMaterializer
import cromiam.cromwell.CromwellClient
import cromiam.server.status.StatusService
import org.scalatest._

import scala.concurrent.{ExecutionContextExecutor, Future}

class CromIamApiServiceSpec extends FlatSpec with Matchers with CromIamApiService with ScalatestRouteTest {
  override def testConfigSource = "akka.loglevel = DEBUG"
  val logger = NoLogging
  // We need an auth object, but it's never used, so...
  val emptyAuthObject = new Authorization(null)

  override def configuration = throw new NotImplementedError("This spec shouldn't need to access the real interface/port")

  override lazy val cromwellClient = new MockCromwellClient()

  override val statusService: StatusService = new StatusService(() => Map.empty)

  "Stats endpoint" should "be forbidden" in {
    Get("/engine/v1/stats").withHeaders(new Authorization(null)) ~> allRoutes ~> check {
      status shouldBe Forbidden
      response.entity shouldBe CromIamApiService.CromIamStatsForbidden.entity
    }
  }
}

class MockCromwellClient()(implicit system: ActorSystem,
                           ece: ExecutionContextExecutor,
                           materializer: ActorMaterializer)
  extends CromwellClient("foo", "bar", 1)(system, ece, materializer) {
  override def forwardToCromwell(httpRequest: HttpRequest): Future[HttpResponse] = {
    Future.successful(HttpResponse(status = InternalServerError))
  }
}
