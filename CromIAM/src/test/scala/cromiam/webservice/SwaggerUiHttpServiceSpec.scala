package cromiam.webservice

import akka.http.scaladsl.model.headers.Location
import akka.http.scaladsl.model.{ContentTypes, StatusCodes, Uri}
import akka.http.scaladsl.server.Route
import akka.http.scaladsl.testkit.ScalatestRouteTest
import common.assertion.CromwellTimeoutSpec
import cromiam.server.config.SwaggerOauthConfig
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

trait SwaggerUiHttpServiceSpec
    extends AnyFlatSpec
    with CromwellTimeoutSpec
    with Matchers
    with ScalatestRouteTest
    with SwaggerUiHttpService {
  def actorRefFactory = system
}

trait SwaggerResourceHttpServiceSpec
    extends AnyFlatSpec
    with CromwellTimeoutSpec
    with Matchers
    with ScalatestRouteTest
    with TableDrivenPropertyChecks
    with SwaggerResourceHttpService {

  val testPathsForOptions = Table("endpoint",
                                  "/",
                                  "/swagger",
                                  "/swagger/index.html",
                                  "/api",
                                  "/api/example",
                                  "/api/example?with=param",
                                  "/api/example/path"
  )
}

trait SwaggerUiResourceHttpServiceSpec
    extends SwaggerUiHttpServiceSpec
    with SwaggerResourceHttpServiceSpec
    with SwaggerUiResourceHttpService

object SwaggerUiHttpServiceSpec {
  // TODO: Re-common-ize swagger out of cromwell's engine and reuse.
  val TestSwaggerUiVersion = "3.23.11" // scala-steward:off
  val SwaggerIndexPreamble =
    """
      |<!-- HTML for static distribution bundle build -->
      |<!DOCTYPE html>
      |<html lang="en">
      |  <head>
      |""".stripMargin.trim
}

class BasicSwaggerUiHttpServiceSpec extends SwaggerUiHttpServiceSpec {
  behavior of "SwaggerUiHttpService"

  override protected def rewriteSwaggerIndex(data: String): String =
    // Replace same magic string used in SwaggerUiResourceHttpService.rewriteSwaggerIndex
    data.replace("window.ui = ui", "replaced-client-id")

  it should "redirect /swagger to /" in {
    Get("/swagger") ~> swaggerUiRoute ~> check {
      status should be(StatusCodes.MovedPermanently)
      header("Location") should be(Option(Location(Uri("/"))))
    }
  }

  it should "redirect /swagger/index.html?url=/swagger/cromiam.yaml to /" in {
    Get("/swagger/index.html?url=/swagger/cromiam.yaml") ~> swaggerUiRoute ~> check {
      status should be(StatusCodes.MovedPermanently)
      header("Location") should be(Option(Location(Uri("/"))))
      responseAs[String] shouldEqual """This and all future requests should be directed to <a href="/">this URI</a>."""
    }
  }

  it should "not return options for /" in {
    Options() ~> Route.seal(swaggerUiRoute) ~> check {
      status should be(StatusCodes.MethodNotAllowed)
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return index.html from the swagger-ui jar for /" in {
    Get("/") ~> swaggerUiRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String].take(15) should be("<!-- HTML for s")
    }
  }

  it should "return index.html from the swagger-ui jar for base url" in {
    Get() ~> swaggerUiRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String].take(15) should be("<!-- HTML for s")
    }
  }

  override def oauthConfig: SwaggerOauthConfig = SwaggerOauthConfig(
    clientId = "test-client-id",
    realm = "test-realm",
    appName = "test-appname"
  )
}

class YamlSwaggerResourceHttpServiceSpec extends SwaggerResourceHttpServiceSpec {
  override def swaggerServiceName = "testservice"

  behavior of "SwaggerResourceHttpService"

  it should "service swagger yaml" in {
    Get("/swagger/testservice.yaml") ~> swaggerResourceRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String] should startWith("swagger: '2.0'\n")
      contentType should be(ContentTypes.`application/octet-stream`)
    }
  }

  it should "not service swagger json" in {
    Get("/swagger/testservice.json") ~> Route.seal(swaggerResourceRoute) ~> check {
      status should be(StatusCodes.NotFound)
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "not service /swagger" in {
    Get("/swagger") ~> Route.seal(swaggerResourceRoute) ~> check {
      status should be(StatusCodes.NotFound)
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return options for all routes" in {
    forAll(testPathsForOptions) { path =>
      Options(path) ~> swaggerResourceRoute ~> check {
        status should be(StatusCodes.OK)
        responseAs[String] should be("OK")
        contentType should be(ContentTypes.`text/plain(UTF-8)`)
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
      contentType should be(ContentTypes.`application/json`)
    }
  }

  it should "not service swagger yaml" in {
    Get("/swagger/testservice.yaml") ~> Route.seal(swaggerResourceRoute) ~> check {
      status should be(StatusCodes.NotFound)
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "not service /swagger" in {
    Get("/swagger") ~> Route.seal(swaggerResourceRoute) ~> check {
      status should be(StatusCodes.NotFound)
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return options for all routes" in {
    forAll(testPathsForOptions) { path =>
      Options(path) ~> swaggerResourceRoute ~> check {
        status should be(StatusCodes.OK)
        responseAs[String] should be("OK")
        contentType should be(ContentTypes.`text/plain(UTF-8)`)
      }
    }
  }
}

class YamlSwaggerUiResourceHttpServiceSpec extends SwaggerUiResourceHttpServiceSpec {

  override def oauthConfig: SwaggerOauthConfig = SwaggerOauthConfig("clientId", "realm", "appName")
  override def swaggerServiceName = "testservice"

  behavior of "SwaggerUiResourceHttpService"

  it should "service swagger yaml" in {
    Get("/swagger/testservice.yaml") ~> swaggerUiResourceRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String] should startWith("swagger: '2.0'\n")
      contentType should be(ContentTypes.`application/octet-stream`)
    }
  }

  it should "not service swagger json" in {
    Get("/swagger/testservice.json") ~> Route.seal(swaggerUiResourceRoute) ~> check {
      status should be(StatusCodes.NotFound)
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return options for all routes" in {
    forAll(testPathsForOptions) { path =>
      Options(path) ~> swaggerUiResourceRoute ~> check {
        status should be(StatusCodes.OK)
        responseAs[String] should be("OK")
        contentType should be(ContentTypes.`text/plain(UTF-8)`)
      }
    }
  }
}

class JsonSwaggerUiResourceHttpServiceSpec extends SwaggerUiResourceHttpServiceSpec {

  override def oauthConfig: SwaggerOauthConfig = SwaggerOauthConfig("clientId", "realm", "appName")

  override def swaggerServiceName = "testservice"
  override def swaggerResourceType = "json"

  behavior of "SwaggerUiResourceHttpService"

  it should "service swagger json" in {
    Get("/swagger/testservice.json") ~> swaggerUiResourceRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String] should startWith("{\n  \"swagger\": \"2.0\",\n")
      contentType should be(ContentTypes.`application/json`)
    }
  }

  it should "not service swagger yaml" in {
    Get("/swagger/testservice.yaml") ~> Route.seal(swaggerUiResourceRoute) ~> check {
      status should be(StatusCodes.NotFound)
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return options for all routes" in {
    forAll(testPathsForOptions) { path =>
      Options(path) ~> swaggerUiResourceRoute ~> check {
        status should be(StatusCodes.OK)
        responseAs[String] should be("OK")
        contentType should be(ContentTypes.`text/plain(UTF-8)`)
      }
    }
  }
}
