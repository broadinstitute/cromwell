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
          |  "errors": [
          |    "first",
          |    "second \"quoted\""
          |  ],
          |  "message": "Bad \"thing\"\nhappened: é",
          |  "status": "fail"
          |}""".stripMargin
    }
    Get() ~> new RuntimeException("solo").errorRequest(StatusCodes.InternalServerError) ~> check {
      responseAs[String] shouldBe "{\n  \"message\": \"solo\",\n  \"status\": \"error\"\n}"
    }
    Get() ~> AggregatedMessageException("none", List.empty).failRequest(StatusCodes.BadRequest) ~> check {
      responseAs[String] shouldBe "{\n  \"errors\": [\n  ],\n  \"message\": \"none\",\n  \"status\": \"fail\"\n}"
    }
  }
}
