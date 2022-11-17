package cromwell.webservice.routes.wes

import akka.http.scaladsl.model.StatusCodes
import cromwell.services.metadata.MetadataService
import cromwell.services.metadata.MetadataService.{WorkflowQueryFailure, WorkflowQuerySuccess}
import cromwell.webservice.routes.wes.WesState.fromStatusString

case class RunListResponse(runs: List[WesRunStatus], next_page_token: String)

object RunListResponse {
  def fromMetadataQueryResponse(response: MetadataService.MetadataQueryResponse): WesResponse = {

      response match {
      case w: WorkflowQuerySuccess =>
        val runs = w.response.results.toList.map(x => WesRunStatus(x.id, fromStatusString(x.status)))
        WesResponseRunList(runs)
      case f: WorkflowQueryFailure => WesErrorResponse(f.reason.getMessage, StatusCodes.BadRequest.intValue)
    }
  }
}
