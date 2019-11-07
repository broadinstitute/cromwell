package cromwell.api.model

import java.time.OffsetDateTime

import org.scalatest.{FlatSpec, Matchers}
import spray.json._
import cromwell.api.model.CromwellQueryResultJsonSupport._

class CromwellQueryResultJsonFormatterSpec extends FlatSpec with Matchers {

  behavior of "CromwellQueryResultJsonFormat"

  val sampleQueryResult = CromwellQueryResults(results = List(
    CromwellQueryResult(
      Option("switcheroo"),
      WorkflowId.fromString("bee51f36-396d-4e22-8a81-33dedff66bf6"),
      Failed,
      Option(OffsetDateTime.parse("2017-07-24T14:44:34.010Z")),
        Option(OffsetDateTime.parse("2017-07-24T14:44:33.227Z")),
      "Archived"
    ),
    CromwellQueryResult(
      Option("switcheroo"),
      WorkflowId.fromString("0071495e-39eb-478e-bc98-8614b986c91e"),
      Succeeded,
        Option(OffsetDateTime.parse("2017-07-24T15:06:45.940Z")),
          Option(OffsetDateTime.parse("2017-07-24T15:04:54.372Z")),
      "Unarchived"
    ),
  ))

  val sampleJson =  """|{
                       |  "results": [
                       |    {
                       |      "name": "switcheroo",
                       |      "id": "bee51f36-396d-4e22-8a81-33dedff66bf6",
                       |      "status": "Failed",
                       |      "end": "2017-07-24T14:44:34.010Z",
                       |      "start": "2017-07-24T14:44:33.227Z",
                       |      "metadataArchiveStatus": "Archived"
                       |    },
                       |    {
                       |      "name": "switcheroo",
                       |      "id": "0071495e-39eb-478e-bc98-8614b986c91e",
                       |      "status": "Succeeded",
                       |      "end": "2017-07-24T15:06:45.940Z",
                       |      "start": "2017-07-24T15:04:54.372Z",
                       |      "metadataArchiveStatus": "Unarchived"
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
