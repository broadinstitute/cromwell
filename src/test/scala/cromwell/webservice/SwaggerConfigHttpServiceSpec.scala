package cromwell.webservice

import com.typesafe.config.ConfigFactory
import org.scalatest.{FlatSpec, Matchers}
import spray.http.StatusCodes
import spray.testkit.ScalatestRouteTest

class SwaggerConfigHttpServiceSpec extends FlatSpec with SwaggerConfigHttpService with ScalatestRouteTest with Matchers {
  override def actorRefFactory = system

  override def apiTypes = Nil

  // Version values in application.conf, possibly reference.conf,
  // and this test conf must be manually synced with build.sbt
  override def swaggerConfig = ConfigFactory.parseString(
    """
      |baseUrl = ""
      |api.version = "1.2.3"
      |ui.version = "2.1.1"
    """.stripMargin)

  behavior of "SwaggerConfigHttpService"

  it should "redirect to the main page when requesting /swagger" in {
    Get("/swagger") ~>
      uiRoutes ~>
      check {
        assertResult(StatusCodes.TemporaryRedirect) {
          status
        }
      }
  }

  it should "return index.html from the swagger-ui jar" in {
    Get("/swagger/index.html") ~>
      uiRoutes ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
      }
  }
}
