package wes2cromwell

import spray.json.{JsonFormat, JsonParser }

case class RunListResponse(runs: List[WesRunStatus], next_page_token: String)

object RunListResponse {
  import RunListResponseJsonSupport._

  def toCromwellQueryResponse(json: String): RunListResponse = {
    val jsonAst = JsonParser(json)
    jsonAst.convertTo[RunListResponse]
  }
}

object RunListResponseJsonSupport {
  import WesResponseJsonSupport._

  implicit val cromwellQueryResponseFormat: JsonFormat[RunListResponse] = jsonFormat2(RunListResponse.apply)
}