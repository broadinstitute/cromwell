package lenthall.spray

import lenthall.spray.WrappedRoute._
import org.scalatest.{FlatSpec, Matchers}
import spray.http.StatusCodes
import spray.routing.HttpService
import spray.testkit.ScalatestRouteTest

class WrappedRouteSpec extends FlatSpec with Matchers with HttpService with ScalatestRouteTest {
  override def actorRefFactory = system

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
    Get("/hello") ~> sealRoute(wrappedRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }

    Post("/hello") ~> sealRoute(wrappedRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }
  }

  it should "not service other routes" in {
    Delete("/hello") ~> sealRoute(wrappedRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }

    Get("/swagger") ~> sealRoute(wrappedRoute) ~> check {
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
    Delete("/hello") ~> sealRoute(legacyWrappedRoute) ~> check {
      status should be(StatusCodes.MethodNotAllowed)
    }

    Get("/swagger") ~> sealRoute(legacyWrappedRoute) ~> check {
      status should be(StatusCodes.NotFound)
    }
  }
}
