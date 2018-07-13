package wes2cromwell

import java.time.ZonedDateTime
import java.time.format.DateTimeFormatter

import akka.actor.{ActorRef, ActorSystem}
import akka.http.scaladsl.Http
import akka.http.scaladsl.model._
import akka.stream.scaladsl.Source
import akka.util.ByteString
import akka.http.scaladsl.model.Multipart.FormData.BodyPart

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}
import akka.http.scaladsl.unmarshalling.Unmarshal
import akka.stream.{ActorMaterializer, Materializer}
import wes2cromwell.WorkflowState._
import spray.json.JsObject

class Transmogriphy(implicit system: ActorSystem, ec: ExecutionContext) {

  val cromwellPath = "http://localhost:8000/api/workflows/v1"

  implicit val materializer: Materializer = ActorMaterializer()

  def getWorkflows(replyTo: ActorRef): Unit = {
    val now = ZonedDateTime.now()
    val oneDayAgoString = now.minusDays(1).format(DateTimeFormatter.ISO_INSTANT).replace(":", "%3A")
    val nowString = now.format(DateTimeFormatter.ISO_INSTANT).replace(":", "%3A")

    val url = cromwellPath + s"/query?start=${oneDayAgoString}&end=${nowString}"
    val request = HttpRequest(method = HttpMethods.GET, uri = url)
    val responseFuture = Http().singleRequest(request)
    responseFuture.onComplete {
      case Success(response) => {
        response.status match {
          case StatusCodes.OK => {
            val bodyDataFuture: Future[String] = Unmarshal(response.entity).to[String]

            bodyDataFuture map { bodyData =>
              val statusList: List[WesResponseStatus] = getWorkflowStatusList(bodyData)
              replyTo ! WesResponseWorkflowList(statusList)
            }
          }

          case StatusCodes.BadRequest => replyTo ! WesResponseError("The request is malformed", response.status.intValue())

          case StatusCodes.InternalServerError => replyTo ! WesResponseError("Cromwell server error", response.status.intValue())

          case _ => replyTo ! WesResponseError("Unexpected response status", response.status.intValue())
        }
      }
      case Failure(_) => replyTo ! WesResponseError("Http error", StatusCodes.InternalServerError.intValue)
    }
  }

  def getWorkflow(replyTo: ActorRef, workflowId: String): Unit = {
    val url: String = cromwellPath + "/" + workflowId + "/metadata"
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

          case StatusCodes.BadRequest => replyTo ! WesResponseError("The request is malformed", response.status.intValue())

          case StatusCodes.InternalServerError => replyTo ! WesResponseError("Cromwell server error", response.status.intValue())

          case _ => replyTo ! WesResponseError("Unexpected response status", response.status.intValue())
        }
      }
      case Failure(_) => replyTo ! WesResponseError("Http error", StatusCodes.InternalServerError.intValue)
    }
  }

  def getWorkflowStatus(replyTo: ActorRef, workflowId: String): Unit = {
    val url: String = cromwellPath + "/" + workflowId + "/status"
    val request = HttpRequest(method = HttpMethods.GET, uri = url)
    val responseFuture = Http().singleRequest(request)
    responseFuture.onComplete {
      case Success(response) => {
        response.status match {
          case StatusCodes.OK => {
            val bodyDataFuture: Future[String] = Unmarshal(response.entity).to[String]

            bodyDataFuture map { bodyData =>
              val cromwellStatusResponse: CromwellStatusResponse = CromwellStatusResponse.toCromwellStatusResponse(bodyData)
              replyTo ! WesResponseStatus(cromwellStatusResponse.id, cromwellToWesStatus(cromwellStatusResponse.status))
            }
          }
          case StatusCodes.BadRequest => replyTo ! WesResponseError("The request is malformed", response.status.intValue())

          case StatusCodes.InternalServerError => replyTo ! WesResponseError("Cromwell server error", response.status.intValue())

          case _ => replyTo ! WesResponseError("Unexpected response status", response.status.intValue())
        }
      }
      case Failure(_) => replyTo ! WesResponseError("Http error", StatusCodes.InternalServerError.intValue)
    }
  }

  def deleteWorkflow(replyTo: ActorRef, workflowId: String): Unit = {
    val url: String = cromwellPath + "/" + workflowId + "/abort"
    val request = HttpRequest(method = HttpMethods.POST, uri = url)
    val responseFuture = Http().singleRequest(request)
    responseFuture.onComplete {
      case Success(response) => {
        response.status match {
          case StatusCodes.OK => {
            val bodyDataFuture: Future[String] = Unmarshal(response.entity).to[String]

            bodyDataFuture map { bodyData =>
              val cromwellStatusResponse: CromwellStatusResponse = CromwellStatusResponse.toCromwellStatusResponse(bodyData)
              replyTo ! WesResponseDeleteWorkflowId(cromwellStatusResponse.id)
            }
          }
          case StatusCodes.NotFound => replyTo ! WesResponseError("The requested workflow was not found", response.status.intValue())

          case StatusCodes.BadRequest => replyTo ! WesResponseError("The request is malformed", response.status.intValue())

          case StatusCodes.InternalServerError => replyTo ! WesResponseError("Cromwell server error", response.status.intValue())

          case _ => replyTo ! WesResponseError("Unexpected response status", response.status.intValue())
        }
      }
      case Failure(_) => replyTo ! WesResponseError("Http error", StatusCodes.InternalServerError.intValue)
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

  def cromwellToWesStatus(cromwellState: String): WorkflowState = {
    cromwellState match {
      case "On Hold" => PAUSED
      case "Submitted" => QUEUED
      case "Running" => RUNNING
      case "Aborting" => CANCELED
      case "Aborted" => CANCELED
      case "Succeeded" => COMPLETE
      case "Failed" => EXECUTOR_ERROR
      case _ => UNKNOWN
    }
  }

  def cromwellMetadataToWesWorkflowLog(json: String): WorkflowLog = {
    val metadata = CromwellMetadata.toCromwellMetadata(json)

    val workflowLogData = WorkflowLogEntry(
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

    WorkflowLog(
      workflow_id = metadata.id,
      state = cromwellToWesStatus(metadata.status),
      workflow_log = Option(workflowLogData),
      task_logs = Option(taskLogs),
      outputs = metadata.outputs
    )
  }

  def cromwellCallsMetadataEntryToLogEntry(taskName: String, callsMetadata: CromwellCallsMetadata): WorkflowLogEntry = {
    val newTaskName = callsMetadata.shardIndex map {
      case -1 => taskName
      case notMinusOne => s"$taskName.$notMinusOne"
    } getOrElse taskName

    WorkflowLogEntry(
      name = Option(newTaskName),
      cmd = None,
      start_time = callsMetadata.start,
      end_time = callsMetadata.end,
      stdout = callsMetadata.stdout,
      stderr = callsMetadata.stderr,
      exit_code = callsMetadata.returnCode
    )
  }

  def getWorkflowStatusList(bodyData: String): List[WesResponseStatus] = {
    val cromwellQueryResponse: CromwellQueryResponse = CromwellQueryResponse.toCromwellQueryResponse(bodyData)
    val cromwellList: List[CromwellStatusResponse] = cromwellQueryResponse.results
    cromwellList.map(x => WesResponseStatus(x.id, cromwellToWesStatus(x.status)))
  }
}


object Transmogriphy {
  // Quick hack to get around compiler flag for porting code from wes2cromwell. Didn't want to spend much time on this
  // as the whole submission scheme is changing as per #3841
  sealed trait WesBodyPart
  case class WesBodyPartError(msg: String, status_code: Int) extends WesBodyPart
  case class ValidWesBodyPart(bodyPart: BodyPart) extends WesBodyPart
}