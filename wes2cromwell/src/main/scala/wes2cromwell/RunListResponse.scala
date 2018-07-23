package wes2cromwell

import spray.json.{ DefaultJsonProtocol, JsonFormat, JsonParser }

case class RunListResponse(runs: List[WesRunStatus], next_page_token: String)

object RunListResponse {
  import DefaultJsonProtocol._
  implicit val cromwellQueryResponseFormat: JsonFormat[RunListResponse] = jsonFormat2(RunListResponse.apply)

  def toCromwellQueryResponse(json: String): RunListResponse = {
    val jsonAst = JsonParser(json)
    jsonAst.convertTo[RunListResponse]
  }
}

