package cromiam.webservice

import akka.event.{LoggingAdapter, NoLogging}
import akka.http.scaladsl.model._
import akka.http.scaladsl.model.headers.{Authorization, OAuth2BearerToken, RawHeader}
import akka.http.scaladsl.testkit.ScalatestRouteTest
import cromwell.util.SampleWdl.HelloWorld
import org.scalatest.{FlatSpec, Matchers}
import spray.json.DefaultJsonProtocol._
import spray.json._

class SubmissionSupportSpec extends FlatSpec with Matchers with ScalatestRouteTest with SubmissionSupport {
  override val cromwellClient = new MockCromwellClient()
  override val samClient = new MockSamClient()
  override val log: LoggingAdapter = NoLogging

  val authorization = Authorization(OAuth2BearerToken("my-token"))
  val badAuthHeaders: List[HttpHeader] = List(authorization, RawHeader("OIDC_CLAIM_user_id", samClient.unauthorizedUserCollectionStr))
  val goodAuthHeaders: List[HttpHeader] = List(authorization, RawHeader("OIDC_CLAIM_user_id", samClient.authorizedUserCollectionStr))

  val submitPath: String = "/api/workflows/v1"

  val workflowSource = Multipart.FormData.BodyPart("workflowSource", HttpEntity(HelloWorld.workflowSource()))
  val workflowInputs = Multipart.FormData.BodyPart("workflowInputs", HttpEntity(MediaTypes.`application/json`, HelloWorld.rawInputs.toJson.toString))
  val formData = Multipart.FormData(workflowSource, workflowInputs).toEntity()


  "Submit query" should "forward the request to Cromwell for authorized SAM user" in {
    Post(submitPath).withHeaders(goodAuthHeaders).withEntity(formData) ~> submitRoute ~> check {
      status shouldEqual StatusCodes.OK
    }
  }

  it should "fail with BadRequest for unauthorized SAM user and not forward the request to Cromwell" in {
    Post(submitPath).withHeaders(badAuthHeaders).withEntity(formData) ~> submitRoute ~> check {
      status shouldEqual StatusCodes.BadRequest
    }
  }
}
