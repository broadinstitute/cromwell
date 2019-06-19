package cromiam.webservice

import akka.event.{LoggingAdapter, NoLogging}
import akka.http.scaladsl.model.headers.{Authorization, OAuth2BearerToken, RawHeader}
import akka.http.scaladsl.model.{ContentTypes, HttpEntity, HttpHeader, StatusCodes}
import akka.http.scaladsl.testkit.ScalatestRouteTest
import cromiam.auth.{Collection, User}
import cromiam.webservice.QuerySupport.LabelContainsOrException
import org.broadinstitute.dsde.workbench.model.WorkbenchUserId
import org.scalatest.{FlatSpec, Matchers}

class QuerySupportSpec extends FlatSpec with Matchers with ScalatestRouteTest with QuerySupport {
  override val cromwellClient = new MockCromwellClient()
  override val samClient = new MockSamClient()
  override val log: LoggingAdapter = NoLogging

  val authorization = Authorization(OAuth2BearerToken("my-token"))
  val authHeaders: List[HttpHeader] = List(
    authorization,
    RawHeader("OIDC_CLAIM_user_id", MockSamClient.AuthorizedUserCollectionStr)
  )

  val queryPath = "/api/workflows/v1/query"
  val getQuery = s"$queryPath?status=Submitted&label=foo:bar&label=foo:baz"
  val badGetQuery = s"$queryPath?status=Submitted&labelor=foo:bar&label=foo:baz"

  val goodPostEntity = HttpEntity(ContentTypes.`application/json`,
    """|[
       |  {
       |    "start": "2015-11-01T00:00:00-04:00"
       |  },
       |  {
       |    "end": "2015-11-04T00:00:00-04:00"
       |  },
       |  {
       |    "status": "Failed"
       |  },
       |  {
       |    "status": "Succeeded"
       |  },
       |  {
       |    "label": "foo:bar"
       |  },
       |  {
       |    "label": "qwe:ert"
       |  }
       |]
       |""".stripMargin
  )
  val badPostEntity = HttpEntity(ContentTypes.`application/json`,
    """|[
       |  {
       |    "start": "2015-11-01T00:00:00-04:00"
       |  },
       |  {
       |    "end": "2015-11-04T00:00:00-04:00"
       |  },
       |  {
       |    "status": "Failed"
       |  },
       |  {
       |    "status": "Succeeded"
       |  },
       |  {
       |    "labelor": "do-not-worry:this-is-expected"
       |  },
       |  {
       |    "label": "qwe:ert"
       |  }
       |]
       |""".stripMargin
  )

  "GET query" should "peacefully forward to Cromwell if nothing is untoward" in {
    Get(getQuery).withHeaders(authHeaders) ~> queryGetRoute ~> check {
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
      responseAs[String] shouldBe "Response from Cromwell"
    }
  }

  it should "reject requests that contain labelOrs" in {
    Get(badGetQuery).withHeaders(authHeaders) ~> queryGetRoute ~> check {
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
      status shouldEqual StatusCodes.InternalServerError
    }
  }

  "POST query" should "peacefully forward to Cromwell is nothing is untoward" in {
    Post(queryPath).withHeaders(authHeaders).withEntity(goodPostEntity) ~> queryPostRoute ~> check {
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
      responseAs[String] shouldBe "Response from Cromwell"
    }
  }

  it should "reject requests that contain labelOrs" in {
    Post(queryPath).withHeaders(authHeaders).withEntity(badPostEntity) ~> queryPostRoute ~> check {
      contentType should be(ContentTypes.`text/plain(UTF-8)`)
      status shouldEqual StatusCodes.InternalServerError
    }
  }

  "userCollectionLabels" should "add in a User's collection when retrieving no collections" in {
    val user = User(WorkbenchUserId("123456780"), authorization)
    val collections = Nil
    userCollectionLabels(user, collections).toList should contain theSameElementsInOrderAs
      List("caas-collection-name:123456780")
  }

  it should "add in a User's collections in normal circumstances" in {
    val user = User(WorkbenchUserId("123456780"), authorization)
    val collections = List(Collection("col1"), Collection("col2"))
    userCollectionLabels(user, collections).toList should contain theSameElementsInOrderAs
      List("caas-collection-name:123456780", "caas-collection-name:col1", "caas-collection-name:col2")
  }

  "ensureNoLabelOrs" should "throw when a query contains a LabelOr" in {
    val user = User(WorkbenchUserId("123456780"), authorization)
    val labelOrs = List("foo:bar", "foo:baz")
    the[LabelContainsOrException] thrownBy {
      ensureNoLabelOrs(user, labelOrs)
    } should have message "User 123456780 submitted a labels query containing an OR which CromIAM is blocking: LABELS CONTAIN 'foo:bar' OR LABELS CONTAIN 'foo:baz'"
  }

  it should "not throw an exception when a query does not contain a LabelOr" in {
    val user = User(WorkbenchUserId("123456780"), authorization)
    val labelOrs = Nil
    ensureNoLabelOrs(user, labelOrs)
  }
}
