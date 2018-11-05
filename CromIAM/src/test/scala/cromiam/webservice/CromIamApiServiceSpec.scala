package cromiam.webservice

import akka.event.NoLogging
import akka.http.scaladsl.model.HttpHeader
import akka.http.scaladsl.model.headers.{OAuth2BearerToken, RawHeader}
import akka.http.scaladsl.model.StatusCodes._
import akka.http.scaladsl.model.headers.Authorization
import akka.http.scaladsl.testkit.ScalatestRouteTest
import com.typesafe.config.Config
import cromiam.server.status.StatusService
import org.scalatest._

class CromIamApiServiceSpec extends FlatSpec with Matchers with CromIamApiService with ScalatestRouteTest {
  override def testConfigSource = "akka.loglevel = DEBUG"
  val log = NoLogging

  override def rootConfig: Config = throw new NotImplementedError("This spec shouldn't need to access the real config")
  override def configuration = throw new NotImplementedError("This spec shouldn't need to access the real interface/port")

  override lazy val cromwellClient = new MockCromwellClient()
  override lazy val samClient = new MockSamClient()
  override val statusService: StatusService = new StatusService(() => Map.empty)

  val version = "v1"

  val authorization = Authorization(OAuth2BearerToken("my-token"))
  val badAuthHeaders: List[HttpHeader] = List(authorization, RawHeader("OIDC_CLAIM_user_id", cromwellClient.unauthorizedUserCollectionStr))
  val goodAuthHeaders: List[HttpHeader] = List(authorization, RawHeader("OIDC_CLAIM_user_id", cromwellClient.authorizedUserCollectionStr))

  behavior of "Status endpoint"
  it should "return 200 for authorized user who has collection associated with workflow" in {
    Get(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/status").withHeaders(goodAuthHeaders) ~> workflowRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe ""
    }
  }

  it should "return 500 for authorized user who doesn't have collection associated with workflow" in {
    Get(s"/api/workflows/$version/${cromwellClient.workflowIdWithoutCollection}/status").withHeaders(goodAuthHeaders) ~> workflowRoutes ~> check {
      status shouldBe InternalServerError
      responseAs[String] shouldBe s"CromIAM unexpected error: java.lang.IllegalArgumentException: Workflow ${cromwellClient.workflowIdWithoutCollection} has no associated collection"
    }
  }

  it should "return SamDenialException for unauthorized user" in {
    Get(s"/api/workflows/$version/${cromwellClient.subworkflowId}/status").withHeaders(badAuthHeaders) ~> workflowRoutes ~> check {
      status shouldBe Forbidden
      responseAs[String] shouldBe "Access Denied"
    }
  }
}
