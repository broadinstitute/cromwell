package cromiam.auth

import akka.http.scaladsl.testkit.ScalatestRouteTest
import org.scalatest.{FlatSpec, Matchers}
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.Route
import Collection._
import akka.http.scaladsl.model.{ContentTypes, HttpEntity, StatusCodes}

class CollectionSpec extends FlatSpec with Matchers with ScalatestRouteTest {
  val labelStringWithoutCollection = "{\"project\":\"foo\", \"group\":\"bar\"}"
  val labelStringWithCollection = "{\"project\":\"foo\", \"caas-collection-name\":\"test\", \"group\":\"bar\"}"

  val goodEntity = HttpEntity(ContentTypes.`application/json`, labelStringWithoutCollection)
  val badEntity = HttpEntity(ContentTypes.`application/json`, labelStringWithCollection)

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
}
