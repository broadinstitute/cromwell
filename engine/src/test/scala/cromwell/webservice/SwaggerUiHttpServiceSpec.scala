package cromwell.webservice

import akka.http.scaladsl.model.{StatusCodes, Uri}
import akka.http.scaladsl.model.headers.Location
import akka.http.scaladsl.server.Route
import akka.http.scaladsl.testkit.ScalatestRouteTest
import com.typesafe.config.ConfigFactory
import cromwell.webservice.SwaggerUiHttpServiceSpec._
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}


trait SwaggerUiHttpServiceSpec extends FlatSpec with Matchers with ScalatestRouteTest with SwaggerUiHttpService {
  override def swaggerUiVersion = TestSwaggerUiVersion
}

trait SwaggerResourceHttpServiceSpec extends FlatSpec with Matchers with ScalatestRouteTest with
TableDrivenPropertyChecks with SwaggerResourceHttpService {
  val testPathsForOptions = Table("endpoint", "/", "/swagger", "/swagger/index.html", "/api", "/api/example",
    "/api/example?with=param", "/api/example/path")
}

trait SwaggerUiResourceHttpServiceSpec extends SwaggerUiHttpServiceSpec with SwaggerResourceHttpServiceSpec with
SwaggerUiResourceHttpService

object SwaggerUiHttpServiceSpec {
  val TestSwaggerUiVersion = "2.1.1"
  val SwaggerIndexPreamble =
    """
      |<!DOCTYPE html>
      |<html>
      |<head>
      |  <meta charset="UTF-8">
      |  <title>Swagger UI</title>""".stripMargin.trim // workaround IDEA's weird formatting of interpolated strings
}

class BasicSwaggerUiHttpServiceSpec extends SwaggerUiHttpServiceSpec {
  behavior of "SwaggerUiHttpService"

  it should "redirect / to /swagger" in {
    Get() ~> swaggerUiRoute ~> check {
      status should be(StatusCodes.TemporaryRedirect)
      header("Location") should be(Option(Location(Uri("/swagger"))))
    }
  }

  it should "not return options for /" in {
    Options() ~> Route.seal(swaggerUiRoute) ~> check {
      status should be(StatusCodes.MethodNotAllowed)
    }
  }

  it should "redirect /swagger to the index.html" in {
    Get("/swagger") ~> swaggerUiRoute ~> check {
      status should be(StatusCodes.TemporaryRedirect)
      header("Location") should be(Option(Location(Uri("/swagger/index.html?url=/api-docs"))))
    }
  }

  it should "return index.html from the swagger-ui jar" in {
    Get("/swagger/index.html") ~> swaggerUiRoute ~> check {
      status should be(StatusCodes.OK)
    }
  }
}

class NoRedirectRootSwaggerUiHttpServiceSpec extends SwaggerUiHttpServiceSpec {
  override def swaggerUiFromRoot = false

  behavior of "SwaggerUiHttpService"

  it should "not redirect / to /swagger" in {
    Get() ~> Route.seal(swaggerUiRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }
  }

  it should "not return options for /" in {
    Options() ~> Route.seal(swaggerUiRoute) ~> check {
      status should be(StatusCodes.MethodNotAllowed)
    }
  }

  it should "redirect /swagger to the index.html" in {
    Get("/swagger") ~> swaggerUiRoute ~> check {
      status should be(StatusCodes.TemporaryRedirect)
      header("Location") should be(Option(Location(Uri("/swagger/index.html?url=/api-docs"))))
    }
  }

  it should "return index.html from the swagger-ui jar" in {
    Get("/swagger/index.html") ~> swaggerUiRoute ~> check {
      status should be(StatusCodes.OK)
    }
  }
}

class DefaultSwaggerUiConfigHttpServiceSpec extends SwaggerUiHttpServiceSpec with SwaggerUiConfigHttpService {
  override def swaggerUiConfig = ConfigFactory.parseString(s"uiVersion = $TestSwaggerUiVersion")

  behavior of "SwaggerUiConfigHttpService"

  it should "redirect /swagger to the index.html" in {
    Get("/swagger") ~> swaggerUiRoute ~> check {
      status should be(StatusCodes.TemporaryRedirect)
      header("Location") should be(Option(Location(Uri("/swagger/index.html?url=/api-docs"))))
    }
  }

  it should "return index.html from the swagger-ui jar" in {
    Get("/swagger/index.html") ~> swaggerUiRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String].take(SwaggerIndexPreamble.length) should be(SwaggerIndexPreamble)
    }
  }
}

class OverriddenSwaggerUiConfigHttpServiceSpec extends SwaggerUiHttpServiceSpec with SwaggerUiConfigHttpService {
  override def swaggerUiConfig = ConfigFactory.parseString(
    s"""
       |baseUrl = /base
       |docsPath = swagger/lenthall.yaml
       |uiPath = ui/path
       |uiVersion = $TestSwaggerUiVersion
     """.stripMargin)

  behavior of "SwaggerUiConfigHttpService"

  it should "redirect /ui/path to the index.html under /base" in {
    Get("/ui/path") ~> swaggerUiRoute ~> check {
      status should be(StatusCodes.TemporaryRedirect)
      header("Location") should be(Option(Location(Uri("/base/ui/path/index.html?url=/base/swagger/lenthall.yaml"))))
    }
  }

  it should "return index.html from the swagger-ui jar" in {
    Get("/ui/path/index.html") ~> swaggerUiRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String].take(SwaggerIndexPreamble.length) should be(SwaggerIndexPreamble)
    }
  }
}

class YamlSwaggerResourceHttpServiceSpec extends SwaggerResourceHttpServiceSpec {
  override def swaggerServiceName = "testservice"

  behavior of "SwaggerResourceHttpService"

  it should "service swagger yaml" in {
    Get("/swagger/testservice.yaml") ~> swaggerResourceRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String] should startWith("swagger: '2.0'\n")
    }
  }

  it should "not service swagger json" in {
    Get("/swagger/testservice.json") ~> Route.seal(swaggerResourceRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }
  }

  it should "not service /swagger" in {
    Get("/swagger") ~> Route.seal(swaggerResourceRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }
  }

  it should "return options for all routes" in {
    forAll(testPathsForOptions) { path =>
      Options(path) ~> swaggerResourceRoute ~> check {
        status should be(StatusCodes.OK)
        responseAs[String] should be("OK")
      }
    }
  }
}

class JsonSwaggerResourceHttpServiceSpec extends SwaggerResourceHttpServiceSpec {
  override def swaggerServiceName = "testservice"

  override def swaggerResourceType = "json"

  behavior of "SwaggerResourceHttpService"

  it should "service swagger json" in {
    Get("/swagger/testservice.json") ~> swaggerResourceRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String] should startWith("{\n  \"swagger\": \"2.0\",\n")
    }
  }

  it should "not service swagger yaml" in {
    Get("/swagger/testservice.yaml") ~> Route.seal(swaggerResourceRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }
  }

  it should "not service /swagger" in {
    Get("/swagger") ~> Route.seal(swaggerResourceRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }
  }

  it should "return options for all routes" in {
    forAll(testPathsForOptions) { path =>
      Options(path) ~> swaggerResourceRoute ~> check {
        status should be(StatusCodes.OK)
        responseAs[String] should be("OK")
      }
    }
  }
}

class NoOptionsSwaggerResourceHttpServiceSpec extends SwaggerResourceHttpServiceSpec {
  override def swaggerServiceName = "testservice"

  override def swaggerAllOptionsOk = false

  behavior of "SwaggerResourceHttpService"

  it should "service swagger yaml" in {
    Get("/swagger/testservice.yaml") ~> swaggerResourceRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String] should startWith("swagger: '2.0'\n")
    }
  }

  it should "not service swagger json" in {
    Get("/swagger/testservice.json") ~> Route.seal(swaggerResourceRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }
  }

  it should "not service /swagger" in {
    Get("/swagger") ~> Route.seal(swaggerResourceRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }
  }

  it should "not return options for all routes" in {
    forAll(testPathsForOptions) { path =>
      Options(path) ~> Route.seal(swaggerResourceRoute) ~> check {
        status should be(StatusCodes.MethodNotAllowed)
      }
    }
  }
}

class YamlSwaggerUiResourceHttpServiceSpec extends SwaggerUiResourceHttpServiceSpec {
  override def swaggerServiceName = "testservice"

  behavior of "SwaggerUiResourceHttpService"

  it should "redirect /swagger to /swagger/index.html with yaml" in {
    Get("/swagger") ~> swaggerUiResourceRoute ~> check {
      status should be(StatusCodes.TemporaryRedirect)
      header("Location") should be(Option(Location(Uri("/swagger/index.html?url=/swagger/testservice.yaml"))))
    }
  }

  it should "service swagger yaml" in {
    Get("/swagger/testservice.yaml") ~> swaggerUiResourceRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String] should startWith("swagger: '2.0'\n")
    }
  }

  it should "not service swagger json" in {
    Get("/swagger/testservice.json") ~> Route.seal(swaggerUiResourceRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }
  }

  it should "return options for all routes" in {
    forAll(testPathsForOptions) { path =>
      Options(path) ~> swaggerUiResourceRoute ~> check {
        status should be(StatusCodes.OK)
        responseAs[String] should be("OK")
      }
    }
  }
}


class JsonSwaggerUiResourceHttpServiceSpec extends SwaggerUiResourceHttpServiceSpec {
  override def swaggerServiceName = "testservice"

  override def swaggerResourceType = "json"

  behavior of "SwaggerUiResourceHttpService"

  it should "redirect /swagger to /swagger/index.html with yaml with json" in {
    Get("/swagger") ~> swaggerUiResourceRoute ~> check {
      status should be(StatusCodes.TemporaryRedirect)
      header("Location") should be(Option(Location(Uri("/swagger/index.html?url=/swagger/testservice.json"))))
    }
  }

  it should "service swagger json" in {
    Get("/swagger/testservice.json") ~> swaggerUiResourceRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String] should startWith("{\n  \"swagger\": \"2.0\",\n")
    }
  }

  it should "not service swagger yaml" in {
    Get("/swagger/testservice.yaml") ~> Route.seal(swaggerUiResourceRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }
  }

  it should "return options for all routes" in {
    forAll(testPathsForOptions) { path =>
      Options(path) ~> swaggerUiResourceRoute ~> check {
        status should be(StatusCodes.OK)
        responseAs[String] should be("OK")
      }
    }
  }
}
