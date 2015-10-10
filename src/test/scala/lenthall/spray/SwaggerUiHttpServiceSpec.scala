package lenthall.spray

import com.typesafe.config.ConfigFactory
import lenthall.spray.SwaggerUiHttpServiceSpec._
import org.scalatest.{FlatSpec, Matchers}
import spray.http.HttpHeaders.Location
import spray.http.{StatusCodes, Uri}
import spray.testkit.ScalatestRouteTest

trait SwaggerUiHttpServiceSpec extends FlatSpec with Matchers with SwaggerUiHttpService with ScalatestRouteTest

object SwaggerUiHttpServiceSpec {
  val TestUiVersion = "2.1.1"
}

class BasicSwaggerUiHttpServiceSpec extends SwaggerUiHttpServiceSpec {

  behavior of "SwaggerUiHttpService"

  override def actorRefFactory = system

  override val swaggerUiInfo = SwaggerUiInfo(TestUiVersion)

  it should "redirect to the index.html when requesting /swagger" in {
    Get("/swagger") ~> swaggerUiRoutes ~> check {
      status should be(StatusCodes.TemporaryRedirect)
      header("Location") should be(Option(Location(Uri("/swagger/index.html?url=/api-docs"))))
    }
  }

  it should "return index.html from the swagger-ui jar" in {
    Get("/swagger/index.html") ~> swaggerUiRoutes ~> check {
      status should be(StatusCodes.OK)
    }
  }
}

class DefaultConfigSwaggerUiHttpServiceSpec extends SwaggerUiHttpServiceSpec with ConfigSwaggerUiHttpService {

  behavior of "ConfigSwaggerUiHttpService"

  override def actorRefFactory = system

  override def swaggerUiConfig = ConfigFactory.parseString(s"uiVersion = $TestUiVersion")

  it should "redirect to the index.html when requesting /swagger" in {
    Get("/swagger") ~> swaggerUiRoutes ~> check {
      status should be(StatusCodes.TemporaryRedirect)
      header("Location") should be(Option(Location(Uri("/swagger/index.html?url=/api-docs"))))
    }
  }

  it should "return index.html from the swagger-ui jar" in {
    Get("/swagger/index.html") ~> swaggerUiRoutes ~> check {
      status should be(StatusCodes.OK)
    }
  }
}

class OverriddenConfigSwaggerUiHttpServiceSpec extends SwaggerUiHttpServiceSpec with ConfigSwaggerUiHttpService {

  behavior of "ConfigSwaggerUiHttpService"

  override def actorRefFactory = system

  override def swaggerUiConfig = ConfigFactory.parseString(
    s"""
       |baseUrl = /base
       |docsPath = swagger/lenthall.yaml
       |uiPath = ui/path
       |uiVersion = $TestUiVersion
     """.stripMargin)

  it should "redirect to the index.html under /base when requesting /ui/path" in {
    Get("/ui/path") ~> swaggerUiRoutes ~> check {
      status should be(StatusCodes.TemporaryRedirect)
      header("Location") should be(Option(Location(Uri("/base/ui/path/index.html?url=/base/swagger/lenthall.yaml"))))
    }
  }

  it should "return index.html from the swagger-ui jar" in {
    Get("/ui/path/index.html") ~> swaggerUiRoutes ~> check {
      status should be(StatusCodes.OK)
    }
  }
}
