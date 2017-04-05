package cromwell.api

import java.net.URL

import akka.http.scaladsl.Http
import akka.actor.ActorSystem
import akka.http.scaladsl.model.{HttpEntity, _}
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

  lazy val submitEndpoint = s"$cromwellUrl/api/workflows/$apiVersion"
  // Everything else is a suffix off the submit endpoint:
  lazy val batchSubmitEndpoint = s"$submitEndpoint/batch"
  private def workflowSpecificEndpoint(workflowId: WorkflowId, endpoint: String) = s"$submitEndpoint/$workflowId/$endpoint"
  def abortEndpoint(workflowId: WorkflowId) = workflowSpecificEndpoint(workflowId, "abort")
  def statusEndpoint(workflowId: WorkflowId) = workflowSpecificEndpoint(workflowId, "status")
  def metadataEndpoint(workflowId: WorkflowId) = workflowSpecificEndpoint(workflowId, "metadata")
  lazy val backendsEndpoint = s"$submitEndpoint/backends"

  import model.CromwellStatusJsonSupport._
  import model.CromwellBackendsJsonSupport._

  private def requestEntityForSubmit(workflowSubmission: WorkflowSubmission) = {
    val sourceBodyParts = Map(
      "wdlSource" -> Option(workflowSubmission.wdl),
      "workflowInputs" -> workflowSubmission.inputsJson,
      "workflowOptions" -> insertSecrets(workflowSubmission.options, workflowSubmission.refreshToken)
    ) collect { case (name, Some(source)) => Multipart.FormData.BodyPart(name, HttpEntity(MediaTypes.`application/json`, ByteString(source))) }

    val zipBodyParts = Map(
      "wdlDependencies" -> workflowSubmission.zippedImports
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
    val submissionSet = workflow.inputsBatch.map(inputs => WorkflowSingleSubmission(workflow.wdl, Option(inputs), workflow.options, workflow.zippedImports, workflow.refreshToken))

    makeRequest[List[CromwellStatus]](HttpRequest(HttpMethods.POST, batchSubmitEndpoint, List.empty[HttpHeader], requestEntity)) map { statuses =>
      val zipped = submissionSet.zip(statuses)
      zipped map { case (submission, status) =>
        SubmittedWorkflow(WorkflowId.fromString(status.id), cromwellUrl, submission)
      }
    }
  }

  def abort(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowStatus] = getRequest[CromwellStatus](abortEndpoint(workflowId)) map WorkflowStatus.apply
  def status(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowStatus] = getRequest[CromwellStatus](statusEndpoint(workflowId)) map WorkflowStatus.apply
  def metadata(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowMetadata] = getRequest[String](metadataEndpoint(workflowId)) map WorkflowMetadata
  def backends(implicit ec: ExecutionContext): Future[CromwellBackends] = getRequest[CromwellBackends](backendsEndpoint)

  /**
    *
    * @tparam A The type of response expected. Must be supported by an implicit unmarshaller from ResponseEntity.
    */
  private def makeRequest[A](request: HttpRequest)(implicit um: Unmarshaller[ResponseEntity, A], ec: ExecutionContext): Future[A] = for {
    response <- Http().singleRequest(request)
    entity <- Future.fromTry(response.toEntity)
    unmarshalled <- entity.to[A]
  } yield unmarshalled

  private def getRequest[A](uri: String)(implicit um: Unmarshaller[ResponseEntity, A], ec: ExecutionContext): Future[A] = makeRequest[A](HttpRequest(uri = uri))

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
}

object CromwellClient {
  final implicit class EnhancedHttpResponse(val response: HttpResponse) extends AnyVal {

    def toEntity: Try[Unmarshal[ResponseEntity]] = response match {
      case HttpResponse(_: StatusCodes.Success, _, entity, _) => Success(Unmarshal(entity))
      case other => Failure(new UnsuccessfulRequestException(other))
    }
  }

  final case class UnsuccessfulRequestException(httpResponse: HttpResponse) extends Exception {
    override def getMessage = httpResponse.toString
  }
}
