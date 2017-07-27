package cromwell.api

import java.net.URL

import akka.http.scaladsl.Http
import akka.actor.ActorSystem
import akka.http.scaladsl.coding.{Deflate, Gzip, NoCoding}
import akka.http.scaladsl.model.{HttpEntity, _}
import akka.http.scaladsl.model.headers.HttpEncodings
import akka.http.scaladsl.unmarshalling.{Unmarshal, Unmarshaller}
import akka.stream.ActorMaterializer
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.util.ByteString
import cromwell.api.model._
import spray.json._

import scala.concurrent.{ExecutionContext, Future}
import cromwell.api.CromwellClient._

import scala.util.{Failure, Success, Try}

class CromwellClient(val cromwellUrl: URL, val apiVersion: String)(implicit actorSystem: ActorSystem, materializer: ActorMaterializer) {

  lazy val engineEndpoint = s"$cromwellUrl/engine/$apiVersion"
  lazy val submitEndpoint = s"$cromwellUrl/api/workflows/$apiVersion"
  // Everything else is a suffix off the submit endpoint:
  lazy val batchSubmitEndpoint = s"$submitEndpoint/batch"
  private def workflowSpecificEndpoint(workflowId: WorkflowId, endpoint: String) = s"$submitEndpoint/$workflowId/$endpoint"
  def abortEndpoint(workflowId: WorkflowId): String = workflowSpecificEndpoint(workflowId, "abort")
  def statusEndpoint(workflowId: WorkflowId): String = workflowSpecificEndpoint(workflowId, "status")
  def metadataEndpoint(workflowId: WorkflowId): String = workflowSpecificEndpoint(workflowId, "metadata")
  def outputsEndpoint(workflowId: WorkflowId): String = workflowSpecificEndpoint(workflowId, "outputs")
  def logsEndpoint(workflowId: WorkflowId): String = workflowSpecificEndpoint(workflowId, "logs")
  def diffEndpoint(workflowA: WorkflowId, callA: String, indexA: ShardIndex, workflowB: WorkflowId, callB: String, indexB: ShardIndex): String = {
    def shardParam(aOrB: String, s: ShardIndex) = s.index.map(i => s"&index$aOrB=$i.toString").getOrElse("")
    s"$submitEndpoint/callcaching/diff?workflowA=$workflowA&callA=$callA&workflowB=$workflowB&callB=$callB${shardParam("A", indexA)}${shardParam("B", indexB)}"
  }
  lazy val backendsEndpoint = s"$submitEndpoint/backends"
  lazy val versionEndpoint = s"$engineEndpoint/version"

  import model.CromwellStatusJsonSupport._
  import model.WorkflowOutputsJsonSupport._
  import model.WorkflowLogsJsonSupport._
  import model.CromwellBackendsJsonSupport._
  import model.CromwellVersionJsonSupport._
  import model.CallCacheDiffJsonSupport._

  private def requestEntityForSubmit(workflowSubmission: WorkflowSubmission) = {
    import cromwell.api.model.LabelsJsonFormatter._

    val sourceBodyParts = Map(
      "workflowSource" -> Option(workflowSubmission.wdl),
      "workflowType" -> workflowSubmission.workflowType,
      "workflowTypeVersion" -> workflowSubmission.workflowTypeVersion,
      "workflowInputs" -> workflowSubmission.inputsJson,
      "workflowOptions" -> insertSecrets(workflowSubmission.options, workflowSubmission.refreshToken),
      "customLabels" -> Option(workflowSubmission.customLabels.toJson.toString)
    ) collect { case (name, Some(source: String)) => Multipart.FormData.BodyPart(name, HttpEntity(MediaTypes.`application/json`, ByteString(source))) }
    val zipBodyParts = Map(
      "workflowDependencies" -> workflowSubmission.zippedImports
    ) collect { case (name, Some(file)) => Multipart.FormData.BodyPart.fromPath(name, MediaTypes.`application/zip`, file.path) }

    val multipartFormData = Multipart.FormData((sourceBodyParts ++ zipBodyParts).toSeq : _*)
    multipartFormData.toEntity()
  }

  def submit(workflow: WorkflowSubmission)(implicit ec: ExecutionContext): Future[SubmittedWorkflow] = {
    val requestEntity = requestEntityForSubmit(workflow)

    makeRequest[CromwellStatus](HttpRequest(HttpMethods.POST, submitEndpoint, List.empty[HttpHeader], requestEntity)) map { status =>
      SubmittedWorkflow(WorkflowId.fromString(status.id), cromwellUrl, workflow)
    }
  }

  def submitBatch(workflow: WorkflowBatchSubmission)(implicit ec: ExecutionContext): Future[List[SubmittedWorkflow]] = {
    import DefaultJsonProtocol._

    val requestEntity = requestEntityForSubmit(workflow)

    // Make a set of submissions that represent the batch (so we can zip with the results later):
    val submissionSet = workflow.inputsBatch.map(inputs => WorkflowSingleSubmission(
      wdl = workflow.wdl,
      workflowType = workflow.workflowType,
      workflowTypeVersion = workflow.workflowTypeVersion,
      inputsJson = Option(inputs),
      options = workflow.options,
      customLabels = workflow.customLabels,
      zippedImports = workflow.zippedImports,
      refreshToken = workflow.refreshToken))

    makeRequest[List[CromwellStatus]](HttpRequest(HttpMethods.POST, batchSubmitEndpoint, List.empty[HttpHeader], requestEntity)) map { statuses =>
      val zipped = submissionSet.zip(statuses)
      zipped map { case (submission, status) =>
        SubmittedWorkflow(WorkflowId.fromString(status.id), cromwellUrl, submission)
      }
    }
  }

  def abort(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowStatus] = simpleRequest[CromwellStatus](uri = abortEndpoint(workflowId), method = HttpMethods.POST) map WorkflowStatus.apply
  def status(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowStatus] = simpleRequest[CromwellStatus](statusEndpoint(workflowId)) map WorkflowStatus.apply
  def metadata(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowMetadata] = simpleRequest[String](metadataEndpoint(workflowId)) map WorkflowMetadata
  def outputs(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowOutputs] = simpleRequest[WorkflowOutputs](outputsEndpoint(workflowId))
  def logs(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowLogs] = simpleRequest[WorkflowLogsStruct](outputsEndpoint(workflowId)) map WorkflowLogs.apply
  def callCacheDiff(workflowA: WorkflowId, callA: String, shardIndexA: ShardIndex, workflowB: WorkflowId, callB: String, shardIndexB: ShardIndex)(implicit ec: ExecutionContext): Future[CallCacheDiff] =
    simpleRequest[CallCacheDiff](diffEndpoint(workflowA, callA, shardIndexA, workflowB, callB, shardIndexB))
  def backends(implicit ec: ExecutionContext): Future[CromwellBackends] = simpleRequest[CromwellBackends](backendsEndpoint)
  def version(implicit ec: ExecutionContext): Future[CromwellVersion] = simpleRequest[CromwellVersion](versionEndpoint)

  private [api] def executeRequest(request: HttpRequest) = Http().singleRequest(request)

  /**
    *
    * @tparam A The type of response expected. Must be supported by an implicit unmarshaller from ResponseEntity.
    */
  private def makeRequest[A](request: HttpRequest)(implicit um: Unmarshaller[ResponseEntity, A], ec: ExecutionContext): Future[A] = for {
    response <- executeRequest(request)
    decoded <- Future.fromTry(decodeResponse(response))
    entity <- Future.fromTry(decoded.toEntity)
    unmarshalled <- unmarshall(response, entity)(um, ec)
  } yield unmarshalled

  private def unmarshall[A](response: HttpResponse, entity: Unmarshal[ResponseEntity])(implicit um: Unmarshaller[ResponseEntity, A], ec: ExecutionContext): Future[A] = {
    import CromwellFailedResponseExceptionJsonSupport._

    if (response.status.isSuccess()) entity.to[A]
    else entity.to[CromwellFailedResponseException] flatMap Future.failed
  }

  private def simpleRequest[A](uri: String, method: HttpMethod = HttpMethods.GET)(implicit um: Unmarshaller[ResponseEntity, A], ec: ExecutionContext): Future[A] = makeRequest[A](HttpRequest(uri = uri, method = method))

  private def insertSecrets(options: Option[String], refreshToken: Option[String]): Option[String] = {
    import DefaultJsonProtocol._
    val tokenKey = "refresh_token"

    def addToken(optionsMap: Map[String, JsValue]): Map[String, JsValue] = {
      refreshToken match {
        case Some(token) if optionsMap.get(tokenKey).isDefined => optionsMap + (tokenKey -> JsString(token))
        case _ => optionsMap
      }
    }

    options map (o => addToken(o.parseJson.asJsObject.convertTo[Map[String, JsValue]]).toJson.toString)
  }

  private val decoders = Map(
    HttpEncodings.gzip -> Gzip,
    HttpEncodings.deflate -> Deflate,
    HttpEncodings.identity -> NoCoding
  )

  private def decodeResponse(response: HttpResponse): Try[HttpResponse] = {
    decoders.get(response.encoding) map { decoder =>
      Try(decoder.decodeMessage(response))
    } getOrElse Failure(UnsuccessfulRequestException(s"No decoder for ${response.encoding}", response))
  }
}

object CromwellClient {
  final implicit class EnhancedHttpResponse(val response: HttpResponse) extends AnyVal {

    def toEntity: Try[Unmarshal[ResponseEntity]] = response match {
      case HttpResponse(_: StatusCodes.Success, _, entity, _) => Success(Unmarshal(entity))
      case HttpResponse(_: StatusCodes.ServerError, _, entity, _) => Success(Unmarshal(entity))
      case other => Failure(UnsuccessfulRequestException("Unmarshalling error", other))
    }
  }

  final case class UnsuccessfulRequestException(message: String, httpResponse: HttpResponse) extends Exception {
    override def getMessage: String = message + ": " + httpResponse.toString
  }
}
