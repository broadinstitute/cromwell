package wes2cromwell

import java.net.URL
import java.time.ZonedDateTime
import java.time.format.DateTimeFormatter

import akka.actor.{ActorRef, ActorSystem}
import akka.http.scaladsl.Http
import akka.http.scaladsl.model._
import akka.stream.scaladsl.Source
import akka.util.ByteString

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}
import akka.http.scaladsl.unmarshalling.Unmarshal
import akka.stream.{ActorMaterializer, Materializer}
import spray.json.JsObject
import wes2cromwell.Wes2CromwellInterface._

class Wes2CromwellInterface(cromwellPath: URL)(implicit system: ActorSystem, ec: ExecutionContext) {

  val cromwellPathOrig = "http://localhost:8000/api/workflows/v1"
  implicit val materializer: Materializer = ActorMaterializer()

  def submit(submission: WesSubmission, origRequest: HttpRequest): Future[WesResponse] = {
    val cromwellRequest = origRequest
      .copy(uri=cromwellPath.toString)
      .withEntity(submission.entity)
    val cromwellResponse = Http().singleRequest(cromwellRequest)
    cromwellResponse.flatMap({ cr =>
      cr.status match {
        case StatusCodes.Created =>
          Unmarshal(cr.entity).to[String].map(s => WesRunId(CromwellStatusResponse.toCromwellStatusResponse(s).id))
        case StatusCodes.BadRequest => Future.successful(BadRequestError)
        case StatusCodes.InternalServerError => Future.successful(InternalServerError)
        case _ => Future.successful(InternalServerError)
      }
    }).recover({case _ => InternalServerError})
  }

  def abort(workflowId: String, origRequest: HttpRequest): Future[WesResponse] = {
    val cromwellUrl = s"$cromwellPath/$workflowId/abort"
    val cromwellRequest = origRequest.copy(uri=cromwellUrl).withMethod(HttpMethods.POST)
    Http().singleRequest(cromwellRequest).flatMap({ cr =>
      cr.status match {
        case StatusCodes.OK =>
          // FIXME: this bit is not only copy pasta but should be clean up-able
          Unmarshal(cr.entity).to[String].map(s => WesRunId(CromwellStatusResponse.toCromwellStatusResponse(s).id))
        case StatusCodes.BadRequest => Future.successful(NotFoundError) // WES doesn't differentiate between not found & malformed like Cromwell does
        case StatusCodes.Unauthorized => Future.successful(UnauthorizedError)
        case StatusCodes.NotFound => Future.successful(NotFoundError)
        case StatusCodes.Forbidden => Future.successful(ForbiddenError)
        case StatusCodes.InternalServerError => Future.successful(InternalServerError)
        case _ => Future.successful(InternalServerError)
      }
    }).recover({case _ => InternalServerError})
  }

  def status(workflowId: String, origRequest: HttpRequest): Future[WesResponse] = {
    val cromwellUrl = s"$cromwellPath/$workflowId/status"
    val cromwellRequest = origRequest.copy(uri = cromwellUrl)
    Http().singleRequest(cromwellRequest).flatMap({ cr =>
      cr.status match {
        case StatusCodes.OK =>
          Unmarshal(cr.entity).to[String].map(s => {
            CromwellStatusResponse.toCromwellStatusResponse(s).toWesRunStatus
          })
        case StatusCodes.BadRequest => Future.successful(BadRequestError)
        case StatusCodes.NotFound => Future.successful(NotFoundError)
        case StatusCodes.InternalServerError => Future.successful(InternalServerError)
        case _ => Future.successful(InternalServerError)
      }
    })
  }

  def getWorkflows(replyTo: ActorRef): Unit = {
    val now = ZonedDateTime.now()
    val oneDayAgoString = now.minusDays(1).format(DateTimeFormatter.ISO_INSTANT).replace(":", "%3A")
    val nowString = now.format(DateTimeFormatter.ISO_INSTANT).replace(":", "%3A")

    val url = cromwellPathOrig + s"/query?start=${oneDayAgoString}&end=${nowString}"
    val request = HttpRequest(method = HttpMethods.GET, uri = url)
    val responseFuture = Http().singleRequest(request)
    responseFuture.onComplete {
      case Success(response) => {
        response.status match {
          case StatusCodes.OK => {
            val bodyDataFuture: Future[String] = Unmarshal(response.entity).to[String]

            bodyDataFuture map { bodyData =>
              val statusList: List[WesRunStatus] = getWorkflowStatusList(bodyData)
              replyTo ! WesResponseRunList(statusList)
            }
          }

          case StatusCodes.BadRequest => replyTo ! WesErrorResponse("The request is malformed", response.status.intValue())

          case StatusCodes.InternalServerError => replyTo ! WesErrorResponse("Cromwell server error", response.status.intValue())

          case _ => replyTo ! WesErrorResponse("Unexpected response status", response.status.intValue())
        }
      }
      case Failure(_) => replyTo ! WesErrorResponse("Http error", StatusCodes.InternalServerError.intValue)
    }
  }

  def getWorkflow(replyTo: ActorRef, workflowId: String): Unit = {
    val url: String = cromwellPathOrig + "/" + workflowId + "/metadata"
    val request = HttpRequest(method = HttpMethods.GET, uri = url)
    val responseFuture = Http().singleRequest(request)
    responseFuture.onComplete {
      case Success(response) => {
        response.status match {
          case StatusCodes.OK => {
            val bodyDataFuture: Future[String] = Unmarshal(response.entity).to[String]

            bodyDataFuture map { bodyData =>
              replyTo ! WesResponseWorkflowMetadata(cromwellMetadataToWesWorkflowLog(bodyData))
            }
          }

          case StatusCodes.BadRequest => replyTo ! WesErrorResponse("The request is malformed", response.status.intValue())

          case StatusCodes.InternalServerError => replyTo ! WesErrorResponse("Cromwell server error", response.status.intValue())

          case _ => replyTo ! WesErrorResponse("Unexpected response status", response.status.intValue())
        }
      }
      case Failure(_) => replyTo ! WesErrorResponse("Http error", StatusCodes.InternalServerError.intValue)
    }
  }

  def makeJsonEntity(content: String): HttpEntity.Default = {
    val bytes = ByteString(content)
    HttpEntity.Default(ContentTypes.`application/json`, bytes.length.toLong, Source.single(bytes))
  }

  def jsonObjectToEntity(content: JsObject): HttpEntity.Default = {
    val bytes = ByteString(content.toString())
    HttpEntity.Default(ContentTypes.`application/json`, bytes.length.toLong, Source.single(bytes))
  }

  def makeTextEntity(content: String): HttpEntity.Default = {
    val bytes = ByteString(content)
    HttpEntity.Default(ContentTypes.`text/plain(UTF-8)`, bytes.length.toLong, Source.single(bytes))
  }

  def cromwellToWesStatus(cromwellState: String): WesRunState = {
    cromwellState match {
      case "On Hold" => WesRunState.PAUSED
      case "Submitted" => WesRunState.QUEUED
      case "Running" => WesRunState.RUNNING
      case "Aborting" => WesRunState.CANCELED
      case "Aborted" => WesRunState.CANCELED
      case "Succeeded" => WesRunState.COMPLETE
      case "Failed" => WesRunState.EXECUTOR_ERROR
      case _ => WesRunState.UNKNOWN
    }
  }

  def cromwellMetadataToWesWorkflowLog(json: String): WesRunLog = {
    val metadata = CromwellMetadata.fromJson(json)

    val workflowLogData = WesLog(
      name = metadata.workflowName,
      cmd = None,
      start_time = metadata.start,
      end_time = metadata.end,
      stdout = None,
      stderr = None,
      exit_code = None
    )

    val taskLogs = for {
      callsArray <- metadata.calls.toSeq
      (taskName, metadataEntries) <- callsArray
      metadataEntry <- metadataEntries
      logEntry = cromwellCallsMetadataEntryToLogEntry(taskName, metadataEntry)
    } yield logEntry

    WesRunLog(
      run_id = metadata.id,
      state = cromwellToWesStatus(metadata.status),
      run_log = Option(workflowLogData),
      task_logs = Option(taskLogs),
      outputs = metadata.outputs
    )
  }

  def cromwellCallsMetadataEntryToLogEntry(taskName: String, callsMetadata: CromwellCallsMetadata): WesLog = {
    val newTaskName = callsMetadata.shardIndex map {
      case -1 => taskName
      case notMinusOne => s"$taskName.$notMinusOne"
    } getOrElse taskName

    WesLog(
      name = Option(newTaskName),
      cmd = None,
      start_time = callsMetadata.start,
      end_time = callsMetadata.end,
      stdout = callsMetadata.stdout,
      stderr = callsMetadata.stderr,
      exit_code = callsMetadata.returnCode
    )
  }

  def getWorkflowStatusList(bodyData: String): List[WesRunStatus] = {
    val cromwellQueryResponse: CromwellQueryResponse = CromwellQueryResponse.toCromwellQueryResponse(bodyData)
    val cromwellList: List[CromwellStatusResponse] = cromwellQueryResponse.results
    cromwellList.map(x => x.toWesRunStatus)
  }
}

object Wes2CromwellInterface {
  val BadRequestError = WesErrorResponse("The request is malformed", StatusCodes.BadRequest.intValue)
  val InternalServerError = WesErrorResponse("Cromwell server error", StatusCodes.InternalServerError.intValue)
  val UnauthorizedError = WesErrorResponse("The request is unauthorized", StatusCodes.Unauthorized.intValue)
  val ForbiddenError = WesErrorResponse("The requester is not authorized to perform this action", StatusCodes.Forbidden.intValue)
  val NotFoundError = WesErrorResponse("The requested workflow run wasn't found", StatusCodes.NotFound.intValue)
}