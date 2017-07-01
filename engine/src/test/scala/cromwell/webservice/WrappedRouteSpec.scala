package cromwell.webservice

import akka.http.scaladsl.model.StatusCodes
import akka.http.scaladsl.testkit.ScalatestRouteTest
import cromwell.webservice.WrappedRoute._
import org.scalatest.{FlatSpec, Matchers}
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.Route

class WrappedRouteSpec extends FlatSpec with Matchers with ScalatestRouteTest {

  def unwrappedRoute = path("hello") {
    get {
      complete("in get")
    } ~ post {
      complete("in post")
    }
  }

  val wrappedRoute = unwrappedRoute.wrapped("api")

  val legacyWrappedRoute = unwrappedRoute.wrapped("api", routeUnwrapped = true)

  behavior of "WrappedRoute"

  it should "service wrapped routes" in {
    Get("/api/hello") ~> wrappedRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String] should be("in get")
    }

    Post("/api/hello") ~> wrappedRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String] should be("in post")
    }
  }

  it should "not service unwrapped routes" in {
    Get("/hello") ~> Route.seal(wrappedRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }

    Post("/hello") ~> Route.seal(wrappedRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }
  }

  it should "not service other routes" in {
    Delete("/hello") ~> Route.seal(wrappedRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }

    Get("/swagger") ~> Route.seal(wrappedRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }
  }

  it should "service wrapped routes when routeUnwrapped is true" in {
    Get("/api/hello") ~> legacyWrappedRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String] should be("in get")
    }

    Post("/api/hello") ~> legacyWrappedRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String] should be("in post")
    }
  }

  it should "service unwrapped routes when routeUnwrapped is true" in {
    Get("/hello") ~> legacyWrappedRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String] should be("in get")
    }

    Post("/hello") ~> legacyWrappedRoute ~> check {
      status should be(StatusCodes.OK)
      responseAs[String] should be("in post")
    }
  }

  it should "not service other routes when routeUnwrapped is true" in {
    Delete("/hello") ~> Route.seal(legacyWrappedRoute) ~> check {
      status should be(StatusCodes.MethodNotAllowed)
    }

    Get("/swagger") ~> Route.seal(legacyWrappedRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }
  }
}
