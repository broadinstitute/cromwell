package cromwell.webservice.routes.wes

import cromwell.services.metadata.MetadataService
import cromwell.services.metadata.MetadataService.{WorkflowQueryFailure, WorkflowQuerySuccess}
import cromwell.webservice.routes.wes.WesState.fromABC


case class RunListResponse(runs: List[WesRunStatus], next_page_token: String)

object RunListResponse {
  def fromMetadataQueryResponse(response: MetadataService.MetadataQueryResponse): RunListResponse = {
    response match {
      case w: WorkflowQuerySuccess =>
        val runs = w.response.results.toList.map(x => WesRunStatus(x.id, fromABC(x.status)))
        RunListResponse(runs, "Not Yet Implemented")
      case _: WorkflowQueryFailure => ???
    }
  }
}