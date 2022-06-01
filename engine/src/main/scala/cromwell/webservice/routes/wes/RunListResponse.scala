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
  def fromMetadataQueryResponse(response: Option[MetadataService.MetadataQueryResponse]): ErrorOr[RunListResponse] = {
    val status_code = StatusCodes.BadRequest.intValue
    response match {
      case Some(w: WorkflowQuerySuccess) =>
        val runs = w.response.results.toList.map(x => WesRunStatus(x.id, fromStatusString(x.status)))
        RunListResponse(runs, "Not Yet Implemented").validNel
      case Some(w: WorkflowQueryFailure) =>
        val runs = w.reason.errorRequest(StatusCodes.BadRequest)
        RunListResponse(runs, "Not Yet Implemented").validNel
      case None => ().invalidNel
    }
  }
}

(x => WesErrorResponse("Failed to fetch query.", StatusCodes.BadRequest.intValue))

//case w: WorkflowQuerySuccess =>
//  val runs = w.response.results.toList.map(x => WesRunStatus(x.id, fromStatusString(x.status)))
//  ErrorOr[RunListResponse](runs, "Not Yet Implemented")]
//  case _: WorkflowQueryFailure => ??? // FIXME: What do we return when we fail to get a Workflow Query?