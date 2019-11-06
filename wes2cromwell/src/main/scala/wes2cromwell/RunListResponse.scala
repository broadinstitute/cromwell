package wes2cromwell

import cromwell.api.model.CromwellQueryResults
import spray.json.JsonParser

case class RunListResponse(runs: List[WesRunStatus], next_page_token: String)

object RunListResponse {
  def fromJson(json: String): RunListResponse = {
    import cromwell.api.model.CromwellQueryResultJsonSupport._

    val jsonAst = JsonParser(json)
    val queryResults = jsonAst.convertTo[CromwellQueryResults]
    val runs = queryResults.results.toList.map(q => WesRunStatus(q.id.toString,  WesState.fromCromwellStatus(q.status.toString)))
    RunListResponse(runs, "Not Yet Implemented") // FIXME: paging is still a known sore spot
  }
}
