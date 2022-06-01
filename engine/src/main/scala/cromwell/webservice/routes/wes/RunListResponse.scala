package cromwell.webservice.routes.wes

import akka.http.scaladsl.model.StatusCodes
import cats.implicits.catsSyntaxValidatedId
import common.validation.ErrorOr.ErrorOr
import cromwell.services.metadata.MetadataService
import cromwell.services.metadata.MetadataService.{WorkflowQueryFailure, WorkflowQuerySuccess}
import cromwell.webservice.routes.wes.WesState.fromStatusString
import cromwell.webservice.WebServiceUtils._

import scala.None


case class RunListResponse(runs: List[WesRunStatus], next_page_token: String)

object RunListResponse {
  def fromMetadataQueryResponse(response: Option[MetadataService.MetadataQueryResponse]): WesResponse = {
    response match {
      case Some(w: WorkflowQuerySuccess) =>
        val runs = w.response.results.toList.map(x => WesRunStatus(x.id, fromStatusString(x.status)))
        WesResponseRunList(runs)
      case Some(_: WorkflowQueryFailure) => WesErrorResponse("Unable to fetch metadata query", StatusCodes.BadRequest.intValue)
      case None => WesErrorResponse("Unable to fetch metadata query", StatusCodes.NotFound.intValue)
    }
  }
}

// (x => WesErrorResponse("Failed to fetch query.", StatusCodes.BadRequest.intValue))

//case w: WorkflowQuerySuccess =>
//  val runs = w.response.results.toList.map(x => WesRunStatus(x.id, fromStatusString(x.status)))
//  ErrorOr[RunListResponse](runs, "Not Yet Implemented")]
//  case _: WorkflowQueryFailure => ??? // FIXME: What do we return when we fail to get a Workflow Query?