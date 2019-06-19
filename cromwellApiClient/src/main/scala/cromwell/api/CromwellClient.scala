package cromwell.api

import java.net.URL

import akka.actor.ActorSystem
import akka.http.scaladsl.Http
import akka.http.scaladsl.coding.{Deflate, Gzip, NoCoding}
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model._
import akka.http.scaladsl.model.headers.{Authorization, HttpCredentials, HttpEncodings}
import akka.http.scaladsl.unmarshalling.{Unmarshal, Unmarshaller}
import akka.stream.ActorMaterializer
import akka.util.ByteString
import cats.effect.IO
import cromwell.api.CromwellClient._
import cromwell.api.model._
import spray.json._

import scala.concurrent.ExecutionContext

class CromwellClient(val cromwellUrl: URL,
                     val apiVersion: String,
                     val defaultCredentials: Option[HttpCredentials]=None)
                    (implicit actorSystem: ActorSystem, materializer: ActorMaterializer) {

  lazy val defaultAuthorization: Option[Authorization] = defaultCredentials.map { Authorization(_) }
  lazy val defaultHeaders: List[HttpHeader] = defaultAuthorization.toList

  lazy val engineEndpoint = s"$cromwellUrl/engine/$apiVersion"
  lazy val submitEndpoint = s"$cromwellUrl/api/workflows/$apiVersion"
  // Everything else is a suffix off the submit endpoint:
  lazy val batchSubmitEndpoint = s"$submitEndpoint/batch"

  def abortEndpoint(workflowId: WorkflowId): Uri = workflowSpecificGetEndpoint(submitEndpoint, workflowId, "abort")
  def statusEndpoint(workflowId: WorkflowId): Uri = workflowSpecificGetEndpoint(submitEndpoint, workflowId, "status")
  def metadataEndpoint(workflowId: WorkflowId, args: Option[Map[String, List[String]]] = None): Uri = workflowSpecificGetEndpoint(submitEndpoint, workflowId, "metadata", args)
  def outputsEndpoint(workflowId: WorkflowId): Uri = workflowSpecificGetEndpoint(submitEndpoint, workflowId, "outputs")
  def labelsEndpoint(workflowId: WorkflowId): Uri = workflowSpecificGetEndpoint(submitEndpoint, workflowId, "labels")
  def logsEndpoint(workflowId: WorkflowId): Uri = workflowSpecificGetEndpoint(submitEndpoint, workflowId, "logs")
  def diffEndpoint(workflowA: WorkflowId, callA: String, indexA: ShardIndex, workflowB: WorkflowId, callB: String, indexB: ShardIndex): String = {
    def shardParam(aOrB: String, s: ShardIndex) = s.index.map(i => s"&index$aOrB=$i.toString").getOrElse("")
    s"$submitEndpoint/callcaching/diff?workflowA=$workflowA&callA=$callA&workflowB=$workflowB&callB=$callB${shardParam("A", indexA)}${shardParam("B", indexB)}"
  }
  lazy val backendsEndpoint = s"$submitEndpoint/backends"
  lazy val versionEndpoint = s"$engineEndpoint/version"

  import model.CallCacheDiffJsonSupport._
  import model.CromwellBackendsJsonSupport._
  import model.CromwellStatusJsonSupport._
  import model.CromwellVersionJsonSupport._
  import model.WorkflowLabelsJsonSupport._
  import model.WorkflowLogsJsonSupport._
  import model.WorkflowOutputsJsonSupport._

  def submit(workflow: WorkflowSubmission)
            (implicit ec: ExecutionContext): FailureResponseOrT[SubmittedWorkflow] = {
    val requestEntity = requestEntityForSubmit(workflow)

    makeRequest[CromwellStatus](HttpRequest(HttpMethods.POST, submitEndpoint, List.empty[HttpHeader], requestEntity)) map { status =>
      SubmittedWorkflow(WorkflowId.fromString(status.id), cromwellUrl, workflow)
    }
  }

  def submitBatch(workflow: WorkflowBatchSubmission)
                 (implicit ec: ExecutionContext): FailureResponseOrT[List[SubmittedWorkflow]] = {
    import DefaultJsonProtocol._

    val requestEntity = requestEntityForSubmit(workflow)

    // Make a set of submissions that represent the batch (so we can zip with the results later):
    val submissionSet = workflow.inputsBatch.map(inputs => WorkflowSingleSubmission(
      workflowSource = workflow.workflowSource,
      workflowUrl = workflow.workflowUrl,
      workflowRoot = workflow.workflowRoot,
      workflowType = workflow.workflowType,
      workflowTypeVersion = workflow.workflowTypeVersion,
      inputsJson = Option(inputs),
      options = workflow.options,
      labels = workflow.labels,
      zippedImports = workflow.zippedImports))

    makeRequest[List[CromwellStatus]](HttpRequest(HttpMethods.POST, batchSubmitEndpoint, List.empty[HttpHeader], requestEntity)) map { statuses =>
      val zipped = submissionSet.zip(statuses)
      zipped map { case (submission, status) =>
        SubmittedWorkflow(WorkflowId.fromString(status.id), cromwellUrl, submission)
      }
    }
  }

  def abort(workflowId: WorkflowId)(implicit ec: ExecutionContext): FailureResponseOrT[WorkflowStatus] = {
    simpleRequest[CromwellStatus](uri = abortEndpoint(workflowId), method = HttpMethods.POST) map WorkflowStatus.apply
  }

  def status(workflowId: WorkflowId)(implicit ec: ExecutionContext): FailureResponseOrT[WorkflowStatus] = {
    simpleRequest[CromwellStatus](statusEndpoint(workflowId)) map WorkflowStatus.apply
  }

  def metadata(workflowId: WorkflowId,
               args: Option[Map[String, List[String]]] = None,
               headers: List[HttpHeader] = defaultHeaders
              )(implicit ec: ExecutionContext): FailureResponseOrT[WorkflowMetadata] = {
    simpleRequest[String](metadataEndpoint(workflowId, args), headers=headers) map WorkflowMetadata
  }

  def outputs(workflowId: WorkflowId)(implicit ec: ExecutionContext): FailureResponseOrT[WorkflowOutputs] = {
    simpleRequest[WorkflowOutputs](outputsEndpoint(workflowId))
  }

  def labels(workflowId: WorkflowId, headers: List[HttpHeader] = defaultHeaders)
            (implicit ec: ExecutionContext): FailureResponseOrT[WorkflowLabels] = {
    simpleRequest[WorkflowLabels](labelsEndpoint(workflowId), headers=headers)
  }

  def logs(workflowId: WorkflowId)(implicit ec: ExecutionContext): FailureResponseOrT[WorkflowLogs] = {
    simpleRequest[WorkflowLogsStruct](outputsEndpoint(workflowId)) map WorkflowLogs.apply
  }

  def callCacheDiff(workflowA: WorkflowId,
                    callA: String,
                    shardIndexA: ShardIndex,
                    workflowB: WorkflowId,
                    callB: String,
                    shardIndexB: ShardIndex
                   )(implicit ec: ExecutionContext): FailureResponseOrT[CallCacheDiff] = {
    simpleRequest[CallCacheDiff](diffEndpoint(workflowA, callA, shardIndexA, workflowB, callB, shardIndexB))
  }

  def backends(implicit ec: ExecutionContext): FailureResponseOrT[CromwellBackends] = {
    simpleRequest[CromwellBackends](backendsEndpoint)
  }

  def version(implicit ec: ExecutionContext): FailureResponseOrT[CromwellVersion] = {
    simpleRequest[CromwellVersion](versionEndpoint)
  }

  private [api] def executeRequest(request: HttpRequest, headers: List[HttpHeader]) = Http().singleRequest(request.withHeaders(headers))

  /**
    *
    * @tparam A The type of response expected. Must be supported by an implicit unmarshaller from ResponseEntity.
    */
  private def makeRequest[A](request: HttpRequest, headers: List[HttpHeader] = defaultHeaders)
                            (implicit um: Unmarshaller[ResponseEntity, A], ec: ExecutionContext):
  FailureResponseOrT[A] = {
    for {
      response <- executeRequest(request, headers).asFailureResponseOrT
      decoded <- FailureResponseOrT.right(decodeResponse(response))
      entity <- FailureResponseOrT.right(decoded.toEntity)
      unmarshalled <- FailureResponseOrT.right(IO.fromFuture(IO(entity.to[A])))
    } yield unmarshalled
  }

  private def simpleRequest[A](uri: Uri,
                               method: HttpMethod = HttpMethods.GET,
                               headers: List[HttpHeader] = defaultHeaders)
                              (implicit um: Unmarshaller[ResponseEntity, A],
                               ec: ExecutionContext): FailureResponseOrT[A] = {
    makeRequest[A](HttpRequest(uri = uri, method = method), headers)
  }

  private val decoders = Map(
    HttpEncodings.gzip -> Gzip,
    HttpEncodings.deflate -> Deflate,
    HttpEncodings.identity -> NoCoding
  )

  private def decodeResponse(response: HttpResponse): IO[HttpResponse] = {
    decoders.get(response.encoding) map { decoder =>
      IO(decoder.decodeMessage(response))
    } getOrElse IO.raiseError(UnsuccessfulRequestException(s"No decoder for ${response.encoding}", response))
  }
}

object CromwellClient {
  final implicit class EnhancedHttpResponse(val response: HttpResponse) extends AnyVal {

    def toEntity: IO[Unmarshal[ResponseEntity]] = response match {
      case HttpResponse(_: StatusCodes.Success, _, entity, _) => IO(Unmarshal(entity))
      case other => IO.raiseError(UnsuccessfulRequestException("Unmarshalling error", other))
    }
  }

  final case class UnsuccessfulRequestException(message: String, httpResponse: HttpResponse) extends Exception(message)

  /**
    * Optionally replace a json value. Returns the original json if:
    * - The jsonOption is None
    * - The key is not found in the json
    * - The valueOption is None
    *
    * @param jsonOption The optional json
    * @param key The key
    * @param valueOption The optional value
    * @return The json with the modified value, or the original json
    */
  def replaceJson(jsonOption: Option[String], key: String, valueOption: Option[String]): Option[String] = {
    val newJsonOption = for {
      value <- valueOption
      json <- jsonOption
      newJson = replaceJson(json, key, value)
    } yield newJson

    newJsonOption orElse jsonOption
  }

  /**
    * Replace a json value. Returns the original json if the key is not found in the json.
    *
    * @param json The json
    * @param key The key
    * @param value The value
    * @return The json with the modified value, or the original json
    */
  def replaceJson(json: String, key: String, value: String): String = {
    import DefaultJsonProtocol._
    val newJsonOption = for {
      _ <- Option(json)
      if json.contains(key)
      map = json.parseJson.asJsObject.convertTo[Map[String, JsValue]]
      if map.contains(key)
      newMap = map.updated(key, JsString(value))
      newJson = newMap.toJson.toString
    } yield newJson

    newJsonOption getOrElse json
  }

  private[api] def requestEntityForSubmit(workflowSubmission: WorkflowSubmission): MessageEntity = {
    import cromwell.api.model.LabelsJsonFormatter._

    val sourceBodyParts = Map(
      "workflowSource" -> workflowSubmission.workflowSource,
      "workflowUrl" -> workflowSubmission.workflowUrl,
      "workflowRoot" -> workflowSubmission.workflowRoot,
      "workflowType" -> workflowSubmission.workflowType,
      "workflowTypeVersion" -> workflowSubmission.workflowTypeVersion,
      "workflowInputs" -> workflowSubmission.inputsJson,
      "workflowOptions" -> workflowSubmission.options,
      "labels" -> workflowSubmission.labels.map(_.toJson.toString)
    ) collect {
      case (name, Some(source: String)) =>
        Multipart.FormData.BodyPart(name, HttpEntity(MediaTypes.`application/json`, ByteString(source)))
    }

    val zipBodyParts = Map(
      "workflowDependencies" -> workflowSubmission.zippedImports
    ) collect {
      case (name, Some(file)) => Multipart.FormData.BodyPart.fromPath(name, MediaTypes.`application/zip`, file.path)
    }

    val multipartFormData = Multipart.FormData((sourceBodyParts ++ zipBodyParts).toSeq : _*)
    multipartFormData.toEntity()
  }

  /**
    * @param args an optional map of HTTP arguments which will be added to the URL
    */
  private [api] def workflowSpecificGetEndpoint(submitEndpoint: String, workflowId: WorkflowId, endpoint: String, args: Option[Map[String, List[String]]] = None) = {
    val url = s"$submitEndpoint/$workflowId/$endpoint"
    val queryBuilder = Uri.Query.newBuilder
    args.getOrElse(Map.empty).foreach({
      case (key, l) => l.foreach(v => queryBuilder.+=(key -> v))
    })
    val queryResult = queryBuilder.result()
    Uri(url).withQuery(queryResult)
  }
}
