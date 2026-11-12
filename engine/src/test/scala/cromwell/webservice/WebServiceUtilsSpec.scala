package cromwell.webservice

import akka.http.scaladsl.model.StatusCodes
import akka.http.scaladsl.testkit.ScalatestRouteTest
import common.exception.AggregatedMessageException
import cromwell.webservice.WebServiceUtils.EnhancedThrowable
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class WebServiceUtilsSpec extends AnyFlatSpec with ScalatestRouteTest with Matchers {

  behavior of "WebServiceUtils"

  it should "errorResponseGolden" in {
    val ex = AggregatedMessageException("Bad \"thing\"\nhappened: é", List("first", "second \"quoted\""))
    Get() ~> ex.failRequest(StatusCodes.BadRequest) ~> check {
      responseAs[String] shouldBe
        """{
          |  "status": "fail",
          |  "message": "Bad \"thing\"\nhappened: é",
          |  "errors": [
          |    "first",
          |    "second \"quoted\""
          |  ]
          |}""".stripMargin
    }
    Get() ~> new RuntimeException("solo").errorRequest(StatusCodes.InternalServerError) ~> check {
      responseAs[String] shouldBe "{\n  \"status\": \"error\",\n  \"message\": \"solo\"\n}"
    }
    Get() ~> AggregatedMessageException("none", List.empty).failRequest(StatusCodes.BadRequest) ~> check {
      responseAs[String] shouldBe "{\n  \"status\": \"fail\",\n  \"message\": \"none\",\n  \"errors\": [\n  ]\n}"
    }
  }
}
