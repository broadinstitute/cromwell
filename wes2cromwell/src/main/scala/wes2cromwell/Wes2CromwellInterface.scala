package wes2cromwell

import java.net.URL
import java.time.ZonedDateTime
import java.time.format.DateTimeFormatter

import akka.actor.{ActorRef, ActorSystem}
import akka.http.scaladsl.Http
import akka.http.scaladsl.model._


import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}
import akka.http.scaladsl.unmarshalling.Unmarshal
import akka.stream.{ActorMaterializer, Materializer}
import wes2cromwell.Wes2CromwellInterface._

class Wes2CromwellInterface(cromwellPath: URL)(implicit system: ActorSystem, ec: ExecutionContext) {
  val cromwellPathOrig = "http://localhost:8000/api/workflows/v1" // FIXME: remove
  implicit val materializer: Materializer = ActorMaterializer()

  def runWorkflow(submission: WesSubmission, headers: List[HttpHeader]): Future[WesResponse] = {
    val cromwellRequest = HttpRequest(
      method = HttpMethods.POST,
      uri = cromwellPath.toString,
      entity = submission.entity,
      headers = headers
    )

    val cromwellResponse = Http().singleRequest(cromwellRequest)
    cromwellResponse.flatMap({ cr =>
      cr.status match {
        case StatusCodes.Created =>
          Unmarshal(cr.entity).to[String].map(s => WesRunId(WesRunStatus.fromJson(s).run_id))
        case StatusCodes.BadRequest => Future.successful(BadRequestError)
        case StatusCodes.InternalServerError => Future.successful(InternalServerError)
        case _ => Future.successful(InternalServerError)
      }
    }).recover({case _ => InternalServerError})
  }

  def cancelRun(workflowId: String, headers: List[HttpHeader]): Future[WesResponse] = {
    val cromwellUrl = s"$cromwellPath/$workflowId/abort"
    val cromwellRequest = HttpRequest(method = HttpMethods.POST, uri = cromwellUrl, headers = headers)
    Http().singleRequest(cromwellRequest).flatMap({ cr =>
      cr.status match {
        case StatusCodes.OK =>
          // FIXME: this bit is not only copy pasta but should be clean up-able
          Unmarshal(cr.entity).to[String].map(s => WesRunId(WesRunStatus.fromJson(s).run_id))
        case StatusCodes.BadRequest => Future.successful(NotFoundError) // WES doesn't differentiate between not found & malformed like Cromwell does
        case StatusCodes.Unauthorized => Future.successful(UnauthorizedError)
        case StatusCodes.NotFound => Future.successful(NotFoundError)
        case StatusCodes.Forbidden => Future.successful(ForbiddenError)
        case StatusCodes.InternalServerError => Future.successful(InternalServerError)
        case _ => Future.successful(InternalServerError)
      }
    }).recover({case _ => InternalServerError})
  }

  def runStatus(workflowId: String, headers: List[HttpHeader]): Future[WesResponse] = {
    val cromwellUrl = s"$cromwellPath/$workflowId/status"
    val cromwellRequest = HttpRequest(method = HttpMethods.GET, uri = cromwellUrl, headers = headers)
    Http().singleRequest(cromwellRequest).flatMap({ cr =>
      cr.status match {
        case StatusCodes.OK =>
          Unmarshal(cr.entity).to[String].map(s => {
            WesRunStatus.fromJson(s)
          })
        case StatusCodes.BadRequest => Future.successful(BadRequestError)
        case StatusCodes.NotFound => Future.successful(NotFoundError)
        case StatusCodes.InternalServerError => Future.successful(InternalServerError)
        case _ => Future.successful(InternalServerError)
      }
    }).recover({case _ => InternalServerError})
  }

  def runLog(workflowId: String, headers: List[HttpHeader]): Future[WesResponse] = {
    val cromwellUrl = s"$cromwellPath/$workflowId/metadata"
    val cromwellRequest = HttpRequest(method = HttpMethods.GET, uri = cromwellUrl, headers = headers)
    Http().singleRequest(cromwellRequest).flatMap({ cr =>
      cr.status match {
        case StatusCodes.OK =>
          Unmarshal(cr.entity).to[String].map(s => {
            WesResponseWorkflowMetadata(WesRunLog.fromJson(s))
          })
        case StatusCodes.BadRequest => Future.successful(BadRequestError)
        case StatusCodes.InternalServerError => Future.successful(InternalServerError)
        case _ => Future.successful(InternalServerError)
      }
    }).recover({case _ => InternalServerError})
  }

  def listRuns(pageSize: Option[Int], pageToken: Option[String], headers: List[HttpHeader]): Future[WesResponse] = {
    // FIXME: to handle - page_size, page_token
    // FIXME: How to handle next_page_token in response?
    val cromwellUrl = s"$cromwellPath/query"
    val cromwellRequest = HttpRequest(method = HttpMethods.GET, uri = cromwellUrl, headers = headers)
    Http().singleRequest(cromwellRequest).flatMap({ cr =>
      cr.status match {
        case StatusCodes.OK =>
          Unmarshal(cr.entity).to[String].map(s => {
            WesResponseRunList(RunListResponse.toCromwellQueryResponse(s).runs)
          })
        case StatusCodes.BadRequest => Future.successful(BadRequestError)
        case StatusCodes.InternalServerError => Future.successful(InternalServerError)
        case _ => Future.successful(InternalServerError)
      }
    }).recover({case _ => InternalServerError})
  }
}

object Wes2CromwellInterface {
  val BadRequestError = WesErrorResponse("The request is malformed", StatusCodes.BadRequest.intValue)
  val InternalServerError = WesErrorResponse("Cromwell server error", StatusCodes.InternalServerError.intValue)
  val UnauthorizedError = WesErrorResponse("The request is unauthorized", StatusCodes.Unauthorized.intValue)
  val ForbiddenError = WesErrorResponse("The requester is not authorized to perform this action", StatusCodes.Forbidden.intValue)
  val NotFoundError = WesErrorResponse("The requested workflow run wasn't found", StatusCodes.NotFound.intValue)
}