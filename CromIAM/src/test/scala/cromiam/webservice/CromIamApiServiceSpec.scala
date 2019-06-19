package cromiam.webservice

import akka.event.NoLogging
import akka.http.scaladsl.model.StatusCodes._
import akka.http.scaladsl.model.headers.{Authorization, OAuth2BearerToken, RawHeader}
import akka.http.scaladsl.model.{ContentTypes, HttpEntity, HttpHeader}
import akka.http.scaladsl.server.MissingHeaderRejection
import akka.http.scaladsl.testkit.ScalatestRouteTest
import com.typesafe.config.Config
import cromiam.server.status.StatusService
import org.scalatest._

class CromIamApiServiceSpec extends FlatSpec with Matchers with CromIamApiService with ScalatestRouteTest {
  override def testConfigSource = "akka.loglevel = DEBUG"

  val log = NoLogging

  override def rootConfig: Config = throw new UnsupportedOperationException("This spec shouldn't need to access the real config")

  override def configuration = throw new UnsupportedOperationException("This spec shouldn't need to access the real interface/port")

  override lazy val cromwellClient = new MockCromwellClient()
  override lazy val cromwellAbortClient = new MockCromwellClient()
  override lazy val samClient = new MockSamClient()
  override val statusService: StatusService = new StatusService(() => Map.empty)

  val version = "v1"

  val authorization = Authorization(OAuth2BearerToken("my-token"))
  val badAuthHeaders: List[HttpHeader] = List(authorization, RawHeader("OIDC_CLAIM_user_id", cromwellClient.unauthorizedUserCollectionStr))
  val goodAuthHeaders: List[HttpHeader] = List(authorization, RawHeader("OIDC_CLAIM_user_id", cromwellClient.authorizedUserCollectionStr))


  behavior of "Status endpoint"
  it should "return 200 for authorized user who has collection associated with root workflow" in {
    Get(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/status").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe "Response from Cromwell"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return 200 for authorized user who has collection associated with subworkflow" in {
    Get(s"/api/workflows/$version/${cromwellClient.subworkflowId}/status").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe "Response from Cromwell"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return 500 for authorized user who doesn't have collection associated with workflow" in {
    Get(s"/api/workflows/$version/${cromwellClient.workflowIdWithoutCollection}/status").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe InternalServerError
      responseAs[String] shouldBe s"CromIAM unexpected error: java.lang.IllegalArgumentException: Workflow ${cromwellClient.workflowIdWithoutCollection} has no associated collection"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return SamDenialException for user who doesn't have view permissions" in {
    Get(s"/api/workflows/$version/${cromwellClient.subworkflowId}/status").withHeaders(badAuthHeaders) ~> allRoutes ~> check {
      status shouldBe Forbidden
      responseAs[String] shouldBe "Access Denied"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "reject request if it doesn't contain OIDC_CLAIM_user_id in header" in {
    Get(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/status") ~> allRoutes ~> check {
      rejection shouldEqual MissingHeaderRejection("OIDC_CLAIM_user_id")
    }
  }


  behavior of "Outputs endpoint"
  it should "return 200 for authorized user who has collection associated with root workflow" in {
    Get(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/outputs").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe "Response from Cromwell"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return 200 for authorized user who has collection associated with subworkflow" in {
    Get(s"/api/workflows/$version/${cromwellClient.subworkflowId}/outputs").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe "Response from Cromwell"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return 500 for authorized user who doesn't have collection associated with workflow" in {
    Get(s"/api/workflows/$version/${cromwellClient.workflowIdWithoutCollection}/outputs").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe InternalServerError
      responseAs[String] shouldBe s"CromIAM unexpected error: java.lang.IllegalArgumentException: Workflow ${cromwellClient.workflowIdWithoutCollection} has no associated collection"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return SamDenialException for user who doesn't have view permissions" in {
    Get(s"/api/workflows/$version/${cromwellClient.subworkflowId}/outputs").withHeaders(badAuthHeaders) ~> allRoutes ~> check {
      status shouldBe Forbidden
      responseAs[String] shouldBe "Access Denied"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "reject request if it doesn't contain OIDC_CLAIM_user_id in header" in {
    Get(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/outputs") ~> allRoutes ~> check {
      rejection shouldEqual MissingHeaderRejection("OIDC_CLAIM_user_id")
    }
  }


  behavior of "Metadata endpoint"
  it should "return 200 for authorized user who has collection associated with root workflow" in {
    Get(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/metadata").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe "Response from Cromwell"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return 200 for authorized user who has collection associated with subworkflow" in {
    Get(s"/api/workflows/$version/${cromwellClient.subworkflowId}/metadata").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe "Response from Cromwell"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return 500 for authorized user who doesn't have collection associated with workflow" in {
    Get(s"/api/workflows/$version/${cromwellClient.workflowIdWithoutCollection}/metadata").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe InternalServerError
      responseAs[String] shouldBe s"CromIAM unexpected error: java.lang.IllegalArgumentException: Workflow ${cromwellClient.workflowIdWithoutCollection} has no associated collection"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return SamDenialException for user who doesn't have view permissions" in {
    Get(s"/api/workflows/$version/${cromwellClient.subworkflowId}/metadata").withHeaders(badAuthHeaders) ~> allRoutes ~> check {
      status shouldBe Forbidden
      responseAs[String] shouldBe "Access Denied"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "reject request if it doesn't contain OIDC_CLAIM_user_id in header" in {
    Get(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/metadata") ~> allRoutes ~> check {
      rejection shouldEqual MissingHeaderRejection("OIDC_CLAIM_user_id")
    }
  }


  behavior of "Logs endpoint"
  it should "return 200 for authorized user who has collection associated with root workflow" in {
    Get(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/logs").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe "Response from Cromwell"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return 200 for authorized user who has collection associated with subworkflow" in {
    Get(s"/api/workflows/$version/${cromwellClient.subworkflowId}/logs").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe "Response from Cromwell"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return 500 for authorized user who doesn't have collection associated with workflow" in {
    Get(s"/api/workflows/$version/${cromwellClient.workflowIdWithoutCollection}/logs").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe InternalServerError
      responseAs[String] shouldBe s"CromIAM unexpected error: java.lang.IllegalArgumentException: Workflow ${cromwellClient.workflowIdWithoutCollection} has no associated collection"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return SamDenialException for user who doesn't have view permissions" in {
    Get(s"/api/workflows/$version/${cromwellClient.subworkflowId}/logs").withHeaders(badAuthHeaders) ~> allRoutes ~> check {
      status shouldBe Forbidden
      responseAs[String] shouldBe "Access Denied"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "reject request if it doesn't contain OIDC_CLAIM_user_id in header" in {
    Get(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/logs") ~> allRoutes ~> check {
      rejection shouldEqual MissingHeaderRejection("OIDC_CLAIM_user_id")
    }
  }


  behavior of "GET Labels endpoint"
  it should "return 200 for authorized user who has collection associated with root workflow" in {
    Get(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/labels").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe "Response from Cromwell"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return 200 for authorized user who has collection associated with subworkflow" in {
    Get(s"/api/workflows/$version/${cromwellClient.subworkflowId}/labels").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe "Response from Cromwell"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return 500 for authorized user who doesn't have collection associated with workflow" in {
    Get(s"/api/workflows/$version/${cromwellClient.workflowIdWithoutCollection}/labels").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe InternalServerError
      responseAs[String] shouldBe s"CromIAM unexpected error: java.lang.IllegalArgumentException: Workflow ${cromwellClient.workflowIdWithoutCollection} has no associated collection"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return SamDenialException for user who doesn't have view permissions" in {
    Get(s"/api/workflows/$version/${cromwellClient.subworkflowId}/labels").withHeaders(badAuthHeaders) ~> allRoutes ~> check {
      status shouldBe Forbidden
      responseAs[String] shouldBe "Access Denied"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "reject request if it doesn't contain OIDC_CLAIM_user_id in header" in {
    Get(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/labels") ~> allRoutes ~> check {
      rejection shouldEqual MissingHeaderRejection("OIDC_CLAIM_user_id")
    }
  }


  behavior of "PATCH Labels endpoint"
  it should "successfully forward request to Cromwell if nothing is untoward" in {
    val labels = """{"key-1":"foo","key-2":"bar"}"""
    val labelEntity = HttpEntity(ContentTypes.`application/json`, labels)
    Patch(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/labels").withHeaders(goodAuthHeaders).withEntity(labelEntity) ~> allRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe "Response from Cromwell"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return 500 and reject request that contains caas-collection-name label" in {
    val labels = """{"key-1":"foo","caas-collection-name":"bar"}"""
    val labelEntity = HttpEntity(ContentTypes.`application/json`, labels)

    Patch(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/labels").withHeaders(goodAuthHeaders).withEntity(labelEntity) ~> allRoutes ~> check {
      status shouldBe InternalServerError
      responseAs[String] shouldBe "Submitted labels contain the key caas-collection-name, which is not allowed\n"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return 500 and reject request that has invalid labels" in {
    val labels = """"key-1":"foo""""
    val labelEntity = HttpEntity(ContentTypes.`application/json`, labels)

    Patch(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/labels").withHeaders(goodAuthHeaders).withEntity(labelEntity) ~> allRoutes ~> check {
      status shouldBe InternalServerError
      responseAs[String] shouldBe "Labels must be a valid JSON object, received: \"key-1\":\"foo\"\n"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return 500 for authorized user who doesn't have collection associated with workflow" in {
    val labels = """{"key-1":"foo","key-2":"bar"}"""
    val labelEntity = HttpEntity(ContentTypes.`application/json`, labels)
    Patch(s"/api/workflows/$version/${cromwellClient.workflowIdWithoutCollection}/labels").withHeaders(goodAuthHeaders).withEntity(labelEntity) ~> allRoutes ~> check {
      status shouldBe InternalServerError
      responseAs[String] shouldBe s"CromIAM unexpected error: java.lang.IllegalArgumentException: Workflow ${cromwellClient.workflowIdWithoutCollection} has no associated collection"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return SamDenialException for user who doesn't have update permissions" in {
    val labels = """{"key-1":"foo","key-2":"bar"}"""
    val labelEntity = HttpEntity(ContentTypes.`application/json`, labels)
    Patch(s"/api/workflows/$version/${cromwellClient.subworkflowId}/labels").withHeaders(badAuthHeaders).withEntity(labelEntity) ~> allRoutes ~> check {
      status shouldBe Forbidden
      responseAs[String] shouldBe "Access Denied"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "reject request if it doesn't contain OIDC_CLAIM_user_id in header" in {
    val labels = """{"key-1":"foo","key-2":"bar"}"""
    val labelEntity = HttpEntity(ContentTypes.`application/json`, labels)
    Patch(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/labels").withEntity(labelEntity) ~> allRoutes ~> check {
      rejection shouldEqual MissingHeaderRejection("OIDC_CLAIM_user_id")
    }
  }


  behavior of "Backends endpoint"
  it should "successfully forward request to Cromwell if auth header is provided" in {
    Get(s"/api/workflows/$version/backends").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe "Response from Cromwell"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "reject request if it doesn't contain OIDC_CLAIM_user_id in header" in {
    Get(s"/api/workflows/$version/backends") ~> allRoutes ~> check {
      rejection shouldEqual MissingHeaderRejection("OIDC_CLAIM_user_id")
    }
  }


  behavior of "ReleaseHold endpoint"
  it should "return 200 for authorized user who has collection associated with root workflow" in {
    Post(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/releaseHold").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe "Response from Cromwell"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return 500 for authorized user who doesn't have collection associated with workflow" in {
    Post(s"/api/workflows/$version/${cromwellClient.workflowIdWithoutCollection}/releaseHold").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe InternalServerError
      responseAs[String] shouldBe s"CromIAM unexpected error: java.lang.IllegalArgumentException: Workflow ${cromwellClient.workflowIdWithoutCollection} has no associated collection"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return SamDenialException for user who doesn't have update permissions" in {
    Post(s"/api/workflows/$version/${cromwellClient.subworkflowId}/releaseHold").withHeaders(badAuthHeaders) ~> allRoutes ~> check {
      status shouldBe Forbidden
      responseAs[String] shouldBe "Access Denied"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "reject request if it doesn't contain OIDC_CLAIM_user_id in header" in {
    Post(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/releaseHold") ~> allRoutes ~> check {
      rejection shouldEqual MissingHeaderRejection("OIDC_CLAIM_user_id")
    }
  }


  behavior of "Abort endpoint"
  it should "return 200 for authorized user who has collection associated with root workflow" in {
    Post(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/abort").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe "Response from Cromwell"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return 500 for authorized user who doesn't have collection associated with workflow" in {
    Post(s"/api/workflows/$version/${cromwellClient.workflowIdWithoutCollection}/abort").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe InternalServerError
      responseAs[String] shouldBe s"CromIAM unexpected error: java.lang.IllegalArgumentException: Workflow ${cromwellClient.workflowIdWithoutCollection} has no associated collection"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return SamDenialException for user who doesn't have abort permissions" in {
    Post(s"/api/workflows/$version/${cromwellClient.subworkflowId}/abort").withHeaders(badAuthHeaders) ~> allRoutes ~> check {
      status shouldBe Forbidden
      responseAs[String] shouldBe "Access Denied"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "reject request if it doesn't contain OIDC_CLAIM_user_id in header" in {
    Post(s"/api/workflows/$version/${cromwellClient.rootWorkflowIdWithCollection}/abort") ~> allRoutes ~> check {
      rejection shouldEqual MissingHeaderRejection("OIDC_CLAIM_user_id")
    }
  }


  behavior of "CallCacheDiff endpoint"
  it should "return 200 for authorized user who has collection associated with both workflows" in {
    val callCacheDiffParams = s"workflowA=${cromwellClient.rootWorkflowIdWithCollection}&callA=helloCall&workflowB=${cromwellClient.anotherRootWorkflowIdWithCollection}&callB=helloCall"
    Get(s"/api/workflows/$version/callcaching/diff?$callCacheDiffParams").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe OK
      responseAs[String] shouldBe "Response from Cromwell"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return BadRequest if request is malformed" in {
    val callCacheDiffParams = s"workflowA=${cromwellClient.rootWorkflowIdWithCollection}&callA=helloCall"
    Get(s"/api/workflows/$version/callcaching/diff?$callCacheDiffParams").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe BadRequest
      responseAs[String] shouldBe "Must supply both workflowA and workflowB to the /callcaching/diff endpoint"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return 500 for authorized user who has doesn't have collection associated with any one workflow" in {
    val callCacheDiffParams = s"workflowA=${cromwellClient.rootWorkflowIdWithCollection}&callA=helloCall&workflowB=${cromwellClient.workflowIdWithoutCollection}&callB=helloCall"
    Get(s"/api/workflows/$version/callcaching/diff?$callCacheDiffParams").withHeaders(goodAuthHeaders) ~> allRoutes ~> check {
      status shouldBe InternalServerError
      responseAs[String] shouldBe s"CromIAM unexpected error: java.lang.IllegalArgumentException: Workflow ${cromwellClient.workflowIdWithoutCollection} has no associated collection"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "return SamDenialException for user who doesn't have read permissions" in {
    val callCacheDiffParams = s"workflowA=${cromwellClient.rootWorkflowIdWithCollection}&callA=helloCall&workflowB=${cromwellClient.anotherRootWorkflowIdWithCollection}&callB=helloCall"
    Get(s"/api/workflows/$version/callcaching/diff?$callCacheDiffParams").withHeaders(badAuthHeaders) ~> allRoutes ~> check {
      status shouldBe Forbidden
      responseAs[String] shouldBe "Access Denied"
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
    }
  }

  it should "reject request if it doesn't contain OIDC_CLAIM_user_id in header" in {
    val callCacheDiffParams = s"workflowA=${cromwellClient.rootWorkflowIdWithCollection}&callA=helloCall&workflowB=${cromwellClient.anotherRootWorkflowIdWithCollection}&callB=helloCall"
    Get(s"/api/workflows/$version/callcaching/diff?$callCacheDiffParams") ~> allRoutes ~> check {
      rejection shouldEqual MissingHeaderRejection("OIDC_CLAIM_user_id")
    }
  }
}
