package wes2cromwell

import spray.json.{ DefaultJsonProtocol, JsonFormat, JsonParser }

case class CromwellQueryResponse(results: List[CromwellStatusResponse], totalResultsCount: Int)

object CromwellQueryResponse {
  import DefaultJsonProtocol._
  implicit val cromwellQueryResponseFormat: JsonFormat[CromwellQueryResponse] = jsonFormat2(CromwellQueryResponse.apply)

  def toCromwellQueryResponse(json: String): CromwellQueryResponse = {
    val jsonAst = JsonParser(json)
    jsonAst.convertTo[CromwellQueryResponse]
  }
}

