package cromiam.webservice

import akka.http.scaladsl.model.{ContentTypes, HttpHeader}
import akka.http.scaladsl.model.StatusCodes._
import akka.http.scaladsl.model.headers.{Authorization, OAuth2BearerToken, RawHeader}
import akka.http.scaladsl.testkit.ScalatestRouteTest
import org.scalatest._


class WomtoolRouteSupportSpec extends FlatSpec with Matchers with WomtoolRouteSupport with ScalatestRouteTest {

  override lazy val cromwellClient = new MockCromwellClient()

  val authorization = Authorization(OAuth2BearerToken("my-token"))
  val goodAuthHeaders: List[HttpHeader] = List(authorization, RawHeader("OIDC_CLAIM_user_id", cromwellClient.authorizedUserCollectionStr))

  behavior of "Womtool endpoint"
  it should "return 200 when we send in a request from an authenticated user" in {
    Post(s"/api/womtool/v1/describe").withHeaders(goodAuthHeaders) ~> womtoolRoutes ~> check {
      responseAs[String] shouldBe "Response from Cromwell"
      status shouldBe OK
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

}
