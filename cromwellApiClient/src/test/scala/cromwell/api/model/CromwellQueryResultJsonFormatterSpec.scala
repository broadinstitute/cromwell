package cromwell.api.model

import java.time.OffsetDateTime

import org.scalatest.{FlatSpec, Matchers}
import spray.json._
import cromwell.api.model.CromwellQueryResultJsonFormatter._

class CromwellQueryResultJsonFormatterSpec extends FlatSpec with Matchers {

  behavior of "CromwellQueryResultJsonFormat"

  val sampleQueryResult = CromwellQueryResults(results = List(
    CromwellQueryResult("switcheroo", WorkflowId.fromString("bee51f36-396d-4e22-8a81-33dedff66bf6"), Failed, OffsetDateTime.parse("2017-07-24T14:44:34.010-04:00"), OffsetDateTime.parse("2017-07-24T14:44:33.227-04:00")),
    CromwellQueryResult("switcheroo", WorkflowId.fromString("0071495e-39eb-478e-bc98-8614b986c91e"), Succeeded, OffsetDateTime.parse("2017-07-24T15:06:45.940-04:00"), OffsetDateTime.parse("2017-07-24T15:04:54.372-04:00"))
  ))

  val sampleJson =  """|{
                       |  "results": [
                       |    {
                       |      "name": "switcheroo",
                       |      "id": "bee51f36-396d-4e22-8a81-33dedff66bf6",
                       |      "status": "Failed",
                       |      "end": "2017-07-24T14:44:34.010-04:00",
                       |      "start": "2017-07-24T14:44:33.227-04:00"
                       |    },
                       |    {
                       |      "name": "switcheroo",
                       |      "id": "0071495e-39eb-478e-bc98-8614b986c91e",
                       |      "status": "Succeeded",
                       |      "end": "2017-07-24T15:06:45.940-04:00",
                       |      "start": "2017-07-24T15:04:54.372-04:00"
                       |    }
                       |  ]
                       |}""".stripMargin.parseJson.asJsObject

  it should "write a query result as a structured JsObject" in {

    sampleQueryResult.toJson shouldEqual sampleJson
  }

  it should "read a query result as a structured JsObject" in {
    sampleJson.convertTo[CromwellQueryResults] shouldBe sampleQueryResult
  }
}
