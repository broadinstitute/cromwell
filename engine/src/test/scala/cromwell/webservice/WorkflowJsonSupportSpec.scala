package cromwell.webservice

import org.scalatest.{FlatSpec, Matchers}
import spray.json._
import WorkflowJsonSupport._

class WorkflowJsonSupportSpec extends FlatSpec with Matchers {

  val sampleSuccessResponse1 = SuccessResponse("good", "msg", Option(JsArray(Vector(JsString("data1"), JsString("data2")))))
  val sampleSuccessResponse2 = SuccessResponse("good", "msg", None)
  val sampleSuccessResponseJson1 =
    """{
      |  "status": "good",
      |  "message": "msg",
      |  "data": [ "data1", "data2" ]
      |}""".stripMargin.parseJson
  val sampleSuccessResponseJson2 =
    """{
      |  "status": "good",
      |  "message": "msg"
      |}""".stripMargin.parseJson

  it should "correctly JSON-ify success messages" in {
    sampleSuccessResponse1.toJson should be(sampleSuccessResponseJson1)
    sampleSuccessResponse2.toJson should be(sampleSuccessResponseJson2)
  }

  it should "correctly parse success response JSONs" in {
    sampleSuccessResponseJson1.convertTo[SuccessResponse] should be(sampleSuccessResponse1)
    sampleSuccessResponseJson2.convertTo[SuccessResponse] should be(sampleSuccessResponse2)
  }

  val sampleFailureResponse = FailureResponse("bad", "msg", Option(Vector("error1", "error2")))
  val sampleFailureResponseJson =
    """{
      |  "status": "bad",
      |  "message": "msg",
      |  "errors": [ "error1", "error2" ]
      |}""".stripMargin.parseJson

  it should "correctly JSON-ify a failure message" in {
    sampleFailureResponse.toJson should be(sampleFailureResponseJson)
  }

  it should "correctly parse failure response JSON" in {
    sampleFailureResponseJson.convertTo[FailureResponse] should be(sampleFailureResponse)
  }

}
