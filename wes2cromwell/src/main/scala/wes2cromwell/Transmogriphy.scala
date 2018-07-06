package wes2cromwell

import java.time.ZonedDateTime
import java.time.format.DateTimeFormatter

import akka.actor.{ActorRef, ActorSystem}
import akka.http.scaladsl.Http
import akka.http.scaladsl.model._
import akka.http.scaladsl.model.Multipart.FormData
import akka.http.scaladsl.model.Multipart.FormData.BodyPart

import scala.concurrent.{Await, ExecutionContext, Future}
import scala.util.{Failure, Success}
import akka.http.scaladsl.unmarshalling.Unmarshal
import akka.stream.scaladsl.Source
import akka.stream.{ActorMaterializer, Materializer}
import akka.util.ByteString
import wes2cromwell.WorkflowState._

import scala.concurrent.duration._

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
            // TODO: yet another case of "bad thing"
            val bodyDataFuture : Future[String] = Unmarshal(response.entity).to[String]
            val bodyData : String = Await.result(bodyDataFuture, 1.second)
            val cromwellQueryResponse: CromwellQueryResponse = CromwellQueryResponse.toCromwellQueryResponse(bodyData)

            val cromwellList : List[CromwellStatusResponse] = cromwellQueryResponse.results
            val statusList : List[WesResponseStatus] = cromwellList.map(x => cromwellToWesStatus(x))
            replyTo ! WesResponseWorkflowList(statusList)
          }

          case StatusCodes.BadRequest =>
            replyTo ! WesResponseError("The request is malformed", response.status.intValue())

          case StatusCodes.InternalServerError =>
            replyTo ! WesResponseError("Cromwell server error", response.status.intValue())

          case _ =>
            replyTo ! WesResponseError("Unexpected response status", response.status.intValue())
        }
      }
      case Failure(_) =>
        replyTo ! WesResponseError("Http error", StatusCodes.InternalServerError.intValue)
    }


  }

  def postWorkflow(replyTo: ActorRef, workflowRequest: WorkflowRequest): Unit = {
    /*
     * See https://docs.google.com/document/d/11_qHPBbEg3Hr4Vs3lh3dvU0dLl1I2zt6rmNxEkW1U1U/edit#
     * for details on the workflow request mapping.
     */

    // Build the list of parts
    // TODO: would it be better to use a mutable list instead of re-creating each time?
    var parts = List(
      BodyPart("workflowType", makeTextEntity(workflowRequest.workflow_type)),
      BodyPart("workflowTypeVersion", makeTextEntity(workflowRequest.workflow_type_version)),
    )

    if (workflowRequest.workflow_descriptor.isDefined) {
      parts = BodyPart("workflowSource", makeJsonEntity(workflowRequest.workflow_descriptor.get)) :: parts
    }

    // Params are optional as are all of the parts in the Cromwell request that are drawn from
    if (workflowRequest.workflow_params.isDefined) {
      val params = WorkflowParams.toWorkflowParams(workflowRequest.workflow_params.get)

      if (params.workflowOptions.isDefined) {
        parts = BodyPart("workflowOptions", makeJsonEntity(params.workflowOptions.get)) :: parts
      }

      if (params.workflowDependencies.isDefined) {
        parts = BodyPart("workflowOptions", makeTextEntity(params.dependenciesZip().get)) :: parts
      }

      for ((input, index) <- params.workflowInputs.zipWithIndex) {
        val suffix: String = if (index == 0) "" else s"_${index + 1}"
        parts = BodyPart("workflowInputs" + suffix, makeTextEntity(input)) :: parts
      }

/* TODO: putting in this text results in a "malformed" error from Cromwell.
 * Not sure how this is supposed to be presented in the multipart form.
 * Leave it out for now.
 */
/*
      val onHold: String = if (params.workflowOnHold.getOrElse(false)) "true" else "false"
      parts = BodyPart("workflowOnHold", makeTextEntity(onHold)) :: parts
*/
    }

    val formData = FormData(Source(parts))
    val request = HttpRequest(method = HttpMethods.POST, uri = cromwellPath, entity = formData.toEntity)
    val responseFuture = Http().singleRequest(request)

    responseFuture.onComplete {
      case Success(response) => {
        response.status match {
          case StatusCodes.Created => {
            // TODO: another "bad thing". The Unmarshall returns a future, but the whole response is already
            // retrieved, so there should be nothing to wait for. I'm guessing there is a better way to this
            implicit val materializer: Materializer = ActorMaterializer()
            val bodyDataFuture : Future[String] = Unmarshal(response.entity).to[String]
            val bodyData : String = Await.result(bodyDataFuture, 1.second)
            val cromwellPostResponse : CromwellStatusResponse = CromwellStatusResponse.toCromwellStatusResponse(bodyData)
            replyTo ! WesResponseCreateWorkflowId(cromwellPostResponse.id)
          }

          case StatusCodes.BadRequest =>
            replyTo ! WesResponseError("The request is malformed", response.status.intValue())

          case StatusCodes.InternalServerError =>
            replyTo ! WesResponseError("Cromwell server error", response.status.intValue())

          case _ =>
            replyTo ! WesResponseError("Unexpected response status", response.status.intValue())
        }
      }
      case Failure(_) =>
        replyTo ! WesResponseError("Http error", StatusCodes.InternalServerError.intValue)
    }
  }

  def makeJsonEntity(content: String): HttpEntity.Default = {
    val bytes = ByteString(content)
    HttpEntity.Default(ContentTypes.`application/json`, bytes.length.toLong, Source.single(bytes))
  }

  def makeTextEntity(content: String): HttpEntity.Default = {
    val bytes = ByteString(content)
    HttpEntity.Default(ContentTypes.`text/plain(UTF-8)`, bytes.length.toLong, Source.single(bytes))
  }

  def getWorkflow(workflowId: String): Option[WorkflowLog] = {
    Some(WorkflowLog(
      workflowId,
      WorkflowRequest(
        Some("Say hello family"),
        Some("params"),
        "WDL",
        "1.2.3",
        Some("tags are keys and values in some format"),
        Some("engine params"),
        Some("url")
      ),
      WorkflowState.EXECUTOR_ERROR,
      WorkflowLogEntry("Workflow Family", Seq("hello world"), "start", "end", "stdout", "stderr", 99),
      Seq(
        WorkflowLogEntry("Joe", Seq("hello Joe"), "start", "end", "stdout", "stderr", 98),
        WorkflowLogEntry("Jone", Seq("hello Jane"), "start", "end", "stdout", "stderr", 97)
      ),
      "outputs"
    ))
  }

  def getWorkflowStatus(replyTo: ActorRef, workflowId: String): Unit = {
    val url: String = cromwellPath + "/" + workflowId + "/status"
    val request = HttpRequest(method = HttpMethods.GET, uri = url)
    val responseFuture = Http().singleRequest(request)
    responseFuture.onComplete {
      case Success(response) => {
        response.status match {
          case StatusCodes.OK => {
            // TODO: another case of "bad thing"
            val bodyDataFuture : Future[String] = Unmarshal(response.entity).to[String]
            val bodyData : String = Await.result(bodyDataFuture, 1.second)
            val cromwellStatusResponse : CromwellStatusResponse = CromwellStatusResponse.toCromwellStatusResponse(bodyData)
            replyTo ! cromwellToWesStatus(cromwellStatusResponse)
          }

          case StatusCodes.BadRequest =>
            replyTo ! WesResponseError("The request is malformed", response.status.intValue())

          case StatusCodes.InternalServerError =>
            replyTo ! WesResponseError("Cromwell server error", response.status.intValue())

          case _ =>
            replyTo ! WesResponseError("Unexpected response status", response.status.intValue())
        }
      }
      case Failure(_) =>
        replyTo ! WesResponseError("Http error", StatusCodes.InternalServerError.intValue)
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
            // TODO: another case of "bad thing"
            val bodyDataFuture : Future[String] = Unmarshal(response.entity).to[String]
            val bodyData : String = Await.result(bodyDataFuture, 1.second)
            val cromwellStatusResponse : CromwellStatusResponse = CromwellStatusResponse.toCromwellStatusResponse(bodyData)
            replyTo ! WesResponseDeleteWorkflowId(cromwellStatusResponse.id)
          }

          case StatusCodes.NotFound =>
            replyTo ! WesResponseError("The requested workflow was not found", response.status.intValue())

          case StatusCodes.BadRequest =>
            replyTo ! WesResponseError("The request is malformed", response.status.intValue())

          case StatusCodes.InternalServerError =>
            replyTo ! WesResponseError("Cromwell server error", response.status.intValue())

          case _ =>
            replyTo ! WesResponseError("Unexpected response status", response.status.intValue())
        }
      }
      case Failure(_) =>
        replyTo ! WesResponseError("Http error", StatusCodes.InternalServerError.intValue)
    }

  }

  def cromwellToWesStatus(cromwellStatusResponse: CromwellStatusResponse) : WesResponseStatus = {
    val workflowState : WorkflowState =
      cromwellStatusResponse.status match {
        case "On Hold" => PAUSED
        case "Submitted" => INITIALIZING
        case "Running" => RUNNING
        case "Aborting" => CANCELED
        case "Aborted" => CANCELED
        case "Succeeded" => COMPLETE
        case "Failed" => EXECUTOR_ERROR
        case _ => UNKNOWN
      }

    WesResponseStatus(cromwellStatusResponse.id, workflowState)
  }

}
