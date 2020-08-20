package wes2cromwell

import java.net.URL

import akka.actor.ActorSystem
import akka.http.scaladsl.Http
import akka.http.scaladsl.model._

import scala.concurrent.{ExecutionContext, Future}
import akka.http.scaladsl.unmarshalling.Unmarshal
import wes2cromwell.Wes2CromwellInterface._

final class Wes2CromwellInterface(cromwellPath: URL)(implicit system: ActorSystem, ec: ExecutionContext) {
  def runWorkflow(submission: WesSubmission, headers: List[HttpHeader]): Future[WesResponse] = {
    // FIXME - Should be able to get away with these fromJsons by implementing the proper marshalling
    // Because this request has the entity, it's not going through the standard forwardToCromwell method
    val cromwellRequest = HttpRequest(method = HttpMethods.POST, uri = cromwellPath.toString, headers = headers, entity=submission.entity)
    handleCromwellResponse(Http().singleRequest(cromwellRequest), (s: String) => WesRunId(WesRunStatus.fromJson(s).run_id))
  }

  def cancelRun(workflowId: String, headers: List[HttpHeader]): Future[WesResponse] = {
    val cromwellUrl = s"$cromwellPath/$workflowId/abort"
    forwardToCromwell(cromwellUrl, headers, HttpMethods.POST, (s: String) => WesRunId(WesRunStatus.fromJson(s).run_id))
  }

  def runStatus(workflowId: String, headers: List[HttpHeader]): Future[WesResponse] = {
    val cromwellUrl = s"$cromwellPath/$workflowId/status"
    forwardToCromwell(cromwellUrl, headers, HttpMethods.GET, (s: String) => WesRunStatus.fromJson(s))
  }

  def runLog(workflowId: String, headers: List[HttpHeader]): Future[WesResponse] = {
    val cromwellUrl = s"$cromwellPath/$workflowId/metadata"
    forwardToCromwell(cromwellUrl, headers, HttpMethods.GET, (s: String) => WesResponseWorkflowMetadata(WesRunLog.fromJson(s)))
  }

  def listRuns(pageSize: Option[Int], pageToken: Option[String], headers: List[HttpHeader]): Future[WesResponse] = {
    // FIXME: to handle - page_size, page_token
    // FIXME: How to handle next_page_token in response?
    val cromwellUrl = s"$cromwellPath/query"
    forwardToCromwell(cromwellUrl, headers, HttpMethods.GET , (s: String) => WesResponseRunList(RunListResponse.fromJson(s).runs))
  }
}

object Wes2CromwellInterface {
  def forwardToCromwell(url: String, headers: List[HttpHeader], method: HttpMethod,
                        f: String => WesResponse)(implicit system: ActorSystem, ec: ExecutionContext): Future[WesResponse] = {
    val cromwellRequest = HttpRequest(method = method, uri = url, headers = headers)
    handleCromwellResponse(Http().singleRequest(cromwellRequest), f)
  }

  def handleCromwellResponse(response: Future[HttpResponse], f: String => WesResponse)(implicit ec: ExecutionContext, actorSystem: ActorSystem): Future[WesResponse] = {
    response.flatMap({ cr =>
      cr.status match {
        /*
          Strictly speaking, this is a larger list than what Cromwell typically returns for most endpoints, however
          leaving it here as things like Unauthorized/Forbidden start showing up a lot more in CromIAM which might
          be underneath these requests instead of OG Cromwell
         */
        case StatusCodes.OK | StatusCodes.Created => Unmarshal(cr.entity).to[String].map(s => f(s))
        case StatusCodes.BadRequest => Future.successful(BadRequestError) // In Cromwell, malformed workflow IDs also generate this and that's not supported in WES, but ....
        case StatusCodes.Unauthorized => Future.successful(UnauthorizedError)
        case StatusCodes.NotFound => Future.successful(NotFoundError)
        case StatusCodes.Forbidden => Future.successful(ForbiddenError)
        case StatusCodes.InternalServerError => Future.successful(InternalServerError)
        case _ => Future.successful(InternalServerError)
      }
    }).recover({case _ => InternalServerError})
  }

  // We'll likely want to live in a world where we're giving more info than this, but that world isn't now
  val BadRequestError = WesErrorResponse("The request is malformed", StatusCodes.BadRequest.intValue)
  val InternalServerError = WesErrorResponse("Cromwell server error", StatusCodes.InternalServerError.intValue)
  val UnauthorizedError = WesErrorResponse("The request is unauthorized", StatusCodes.Unauthorized.intValue)
  val ForbiddenError = WesErrorResponse("The requester is not authorized to perform this action", StatusCodes.Forbidden.intValue)
  val NotFoundError = WesErrorResponse("The requested workflow run wasn't found", StatusCodes.NotFound.intValue)
}
