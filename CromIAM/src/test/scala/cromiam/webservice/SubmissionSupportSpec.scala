package cromiam.webservice

import akka.event.{LoggingAdapter, NoLogging}
import akka.http.scaladsl.model._
import akka.http.scaladsl.model.headers.{Authorization, OAuth2BearerToken, RawHeader}
import akka.http.scaladsl.server.AuthorizationFailedRejection
import akka.http.scaladsl.testkit.ScalatestRouteTest
import org.scalatest.{FlatSpec, Matchers}

class SubmissionSupportSpec extends FlatSpec with Matchers with ScalatestRouteTest with SubmissionSupport {
  override val cromwellClient = new MockCromwellClient()
  override val samClient = new MockSamClient()
  override val log: LoggingAdapter = NoLogging

  val authorization = Authorization(OAuth2BearerToken("my-token"))
  val badAuthHeaders: List[HttpHeader] = List(authorization, RawHeader("OIDC_CLAIM_user_id", samClient.unauthorizedUserCollectionStr))
  val goodAuthHeaders: List[HttpHeader] = List(authorization, RawHeader("OIDC_CLAIM_user_id", samClient.authorizedUserCollectionStr))
  val notWhitelistedUserHeader: List[HttpHeader] = List(authorization, RawHeader("OIDC_CLAIM_user_id", samClient.notWhitelistedUser))

  val submitPath: String = "/api/workflows/v1"

  val helloWorldWdl =
    s"""
       |task hello {
       |  String addressee
       |  command {
       |    echo "Hello World!"
       |  }
       |  output {
       |    String salutation = read_string(stdout())
       |  }
       |  RUNTIME
       |}
       |
       |workflow wf_hello {
       |  call hello
       |}
      """.stripMargin
  val workflowSource = Multipart.FormData.BodyPart("workflowSource", HttpEntity(helloWorldWdl))
  val formData = Multipart.FormData(workflowSource).toEntity()


  "Submit endpoint" should "forward the request to Cromwell for authorized SAM user" in {
    Post(submitPath).withHeaders(goodAuthHeaders).withEntity(formData) ~> submitRoute ~> check {
      status shouldEqual StatusCodes.OK
      responseAs[String] shouldBe "Response from Cromwell"
    }
  }

  it should "fail with BadRequest for unauthorized SAM user and not forward the request to Cromwell" in {
    Post(submitPath).withHeaders(badAuthHeaders).withEntity(formData) ~> submitRoute ~> check {
      status shouldEqual StatusCodes.BadRequest
    }
  }

  it should "reject the request when not whitelisted user tries to submit to CromIAM" in {
    Post(submitPath).withHeaders(notWhitelistedUserHeader).withEntity(formData) ~> submitRoute ~> check {
      rejection shouldEqual AuthorizationFailedRejection
    }
  }
}
