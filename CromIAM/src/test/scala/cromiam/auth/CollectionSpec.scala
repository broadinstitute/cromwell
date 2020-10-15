package cromiam.auth

import akka.http.scaladsl.model.{ContentTypes, HttpEntity, StatusCodes}
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.Route
import akka.http.scaladsl.testkit.ScalatestRouteTest
import common.assertion.CromwellTimeoutSpec
import cromiam.auth.Collection._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class CollectionSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with ScalatestRouteTest {
  val labelStringWithoutCollection = "{\"project\":\"foo\", \"group\":\"bar\"}"
  val labelStringWithCollection = "{\"project\":\"foo\", \"caas-collection-name\":\"test\", \"group\":\"bar\"}"

  val goodEntity = HttpEntity(ContentTypes.`application/json`, labelStringWithoutCollection)
  val badEntity = HttpEntity(ContentTypes.`application/json`, labelStringWithCollection)

  val testPath = "test"

  val testRoute: Route = path(testPath) {
    post {
      entity(as[String]) { labelJson =>
        validateLabels(Option(labelJson)) { labels =>
          complete {
            labels.toString
          }
        }
      }
    }
  }

  "The label validation" should "reject requests with a collection label" in {
    Post(s"/$testPath").withEntity(badEntity) ~> testRoute ~> check {
      status shouldBe StatusCodes.InternalServerError
    }
  }

  it should "pass when passed labels which do not include the collection label" in {
    Post(s"/$testPath").withEntity(goodEntity) ~> testRoute ~> check {
      status shouldBe StatusCodes.OK
    }
  }


  behavior of "validateLabels"
  it should "not throw exception when labels are valid" in {
    val labels = """{"key-1":"foo","key-2":"bar"}"""

    validateLabels(Option(labels))
  }

  it should "throw LabelContainsCollectionException when patch labels contain caas-collection-name label" in {
    val labels = """{"key-1":"foo","caas-collection-name":"bar"}"""

    the[LabelContainsCollectionException] thrownBy {
      validateLabels(Option(labels))
    } should have message "Submitted labels contain the key caas-collection-name, which is not allowed\n"
  }

  it should "throw invalid labels exception when labels are invalid" in {
    val labels = """"key-1":"foo""""

    the[InvalidLabelsException] thrownBy {
      validateLabels(Option(labels))
    } should have message "Labels must be a valid JSON object, received: \"key-1\":\"foo\"\n"
  }
}
