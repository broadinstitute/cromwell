import akka.event.NoLogging
import akka.http.scaladsl.model.StatusCodes._
import akka.http.scaladsl.testkit.ScalatestRouteTest
import org.scalatest._
import service.CromIamApiService


class CromIamApiServiceSpec extends FlatSpec with Matchers with ScalatestRouteTest with CromIamApiService {
  override def testConfigSource = "akka.loglevel = WARNING"
  def config = testConfig
  val logger = NoLogging

  "Stats endpoint" should "respond with internal server error" in {
    Get("/api/engine/v1/stats") ~> cromIamRoutes ~> check {
      status shouldBe InternalServerError
    }
  }
}
