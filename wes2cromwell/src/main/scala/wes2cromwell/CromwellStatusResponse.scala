package wes2cromwell

import spray.json.{ DefaultJsonProtocol, JsonFormat, JsonParser }

case class CromwellStatusResponse(id: String, status: String) {
  def toWesRunStatus: WesRunStatus = WesRunStatus(id, NewWesRunState.fromCromwellState(status))
}

object CromwellStatusResponse {
  import DefaultJsonProtocol._
  implicit val cromwellPostResponseFormat: JsonFormat[CromwellStatusResponse] = jsonFormat2(CromwellStatusResponse.apply)

  def toCromwellStatusResponse(json: String): CromwellStatusResponse = {
    val jsonAst = JsonParser(json)
    jsonAst.convertTo[CromwellStatusResponse]
  }
}
