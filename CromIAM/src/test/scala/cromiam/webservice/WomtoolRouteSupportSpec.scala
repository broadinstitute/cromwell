package cromiam.webservice

import akka.http.scaladsl.model.ContentTypes
import akka.http.scaladsl.model.StatusCodes._
import akka.http.scaladsl.model.headers.{Authorization, OAuth2BearerToken, RawHeader}
import akka.http.scaladsl.server.Route.seal
import akka.http.scaladsl.server.{AuthorizationFailedRejection, MissingHeaderRejection}
import akka.http.scaladsl.testkit.ScalatestRouteTest
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class WomtoolRouteSupportSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with WomtoolRouteSupport with ScalatestRouteTest {

  override lazy val cromwellClient = new MockCromwellClient()
  override lazy val samClient = new MockSamClient()

  behavior of "Womtool endpoint routes"

  it should "return 200 when we request to the right path" in {
    Post(
      s"/api/womtool/v1/describe"
    ).withHeaders(
      List(Authorization(OAuth2BearerToken("my-token")), RawHeader("OIDC_CLAIM_user_id", "enabled@example.com"))
    ) ~> womtoolRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe "Hey there, workflow describer"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return 403 when we request with a disabled user" in {
    Post(
      s"/api/womtool/v1/describe"
    ).withHeaders(
      List(Authorization(OAuth2BearerToken("my-token")), RawHeader("OIDC_CLAIM_user_id", "disabled@example.com"))
    ) ~> womtoolRoutes ~> check {
      rejection shouldEqual AuthorizationFailedRejection
    }
  }

  it should "bail out with no user ID" in {
    Post(
      s"/api/womtool/v1/describe"
    ).withHeaders(
      List(Authorization(OAuth2BearerToken("my-token")))
    ) ~> womtoolRoutes ~> check {
      rejection shouldEqual MissingHeaderRejection("OIDC_CLAIM_user_id")
    }
  }

  it should "return 404 when no auth token provided" in {
    Post(
      s"/api/womtool/v1/describe"
    ).withHeaders(
      List(RawHeader("OIDC_CLAIM_user_id", "enabled@example.com"))
      // "[An] explicit call on the Route.seal method is needed in test code, but in your application code it is not necessary."
      // https://doc.akka.io/docs/akka-http/current/routing-dsl/testkit.html#testing-sealed-routes
      // https://doc.akka.io/docs/akka-http/current/routing-dsl/routes.html#sealing-a-route
    ) ~> seal(womtoolRoutes) ~> check {
      responseAs[String] shouldEqual "The requested resource could not be found."
      status shouldBe NotFound
    }
  }

  it should "bail out with no headers" in {
    Post(
      s"/api/womtool/v1/describe"
    ).withHeaders(
      List.empty
    ) ~> womtoolRoutes ~> check {
      rejection shouldEqual MissingHeaderRejection("OIDC_CLAIM_user_id")
    }
  }

}
