package cromiam.webservice

import akka.event.NoLogging
import akka.http.scaladsl.model.{ContentTypes, HttpEntity, HttpHeader, StatusCodes}
import akka.http.scaladsl.model.headers.{Authorization, OAuth2BearerToken, RawHeader}
import akka.http.scaladsl.testkit.ScalatestRouteTest
import cromiam.auth.Collection.LabelContainsCollectionException
import cromiam.auth.{Collection, User}
import org.broadinstitute.dsde.workbench.model.WorkbenchUserId
import org.scalatest.{FlatSpec, Matchers}


class QuerySupportSpec extends FlatSpec with Matchers with ScalatestRouteTest with QuerySupport {
  override val cromwellClient = new MockCromwellClient()
  override val samClient = new MockSamClient()
  override val log = NoLogging

  val authorization = Authorization(OAuth2BearerToken("my-token"))
  val authHeaders: List[HttpHeader] = List(authorization, RawHeader("OIDC_CLAIM_user_id", "123456789"))

  val queryPath = "/api/workflows/v1/query"
  val getQuery = s"$queryPath?status=Submitted&label=foo:bar&label=foo:baz"
  val badGetQuery = s"$queryPath?status=Submitted&label=caas-collection-name:badnews&label=foo:bar&label=foo:baz"

  val goodPostEntity = HttpEntity(ContentTypes.`application/json`, "[{\"start\": \"2015-11-01T00:00:00-04:00\"}, {\"end\": \"2015-11-04T00:00:00-04:00\"}, {\"status\": \"Failed\"}, {\"status\": \"Succeeded\"}, {\"label\": \"foo:bar\"}, {\"label\": \"qwe:ert\"}]")
  val badPostEntity = HttpEntity(ContentTypes.`application/json`, "[{{\"label\": \"caas-collection-name:oops\"}, \"start\": \"2015-11-01T00:00:00-04:00\"}, {\"end\": \"2015-11-04T00:00:00-04:00\"}, {\"status\": \"Failed\"}, {\"status\": \"Succeeded\"}, {\"label\": \"foo:bar\"}, {\"label\": \"qwe:ert\"}]")

  "GET query" should "peacefully forward to Cromwell if nothing is untoward" in {
    Get(getQuery).withHeaders(authHeaders) ~> queryGetRoute ~> check {
      responseAs[String] shouldEqual ""
    }
  }

  it should "reject requests with collection labels specified" in {
    Get(badGetQuery).withHeaders(authHeaders) ~> queryGetRoute ~> check {
      status shouldEqual StatusCodes.InternalServerError
    }
  }

  "POST query" should "peacefully forward to Cromwell is nothing is untoward" in {
    Post(queryPath).withHeaders(authHeaders).withEntity(goodPostEntity) ~> queryPostRoute ~> check {
      responseAs[String] shouldEqual ""
    }
  }

  it should "reject requests with collection labels specified" in {
    Post(queryPath).withHeaders(authHeaders).withEntity(badPostEntity) ~> queryPostRoute ~> check {
      status shouldEqual StatusCodes.InternalServerError
    }
  }

  "processLabels" should "add in a User's collections in normal circumstances" in {
    processLabels(User(WorkbenchUserId("123456780"), authorization),
      List(Collection("col1"), Collection("col2")),
      List("foo:bar", "foo:baz")) shouldBe List("foo:bar", "foo:baz", "caas-collection-name:col1", "caas-collection-name:col2")
  }

  it should "throw when a collection already exists" in {
    assertThrows[LabelContainsCollectionException] {
      processLabels(User(WorkbenchUserId("123456780"), authorization),
        List(Collection("col1"), Collection("col2")),
        List("foo:bar", "foo:baz", "caas-collection-name:oops"))
    }
  }
}
