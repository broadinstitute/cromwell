package centaur.api

import java.util.UUID
import java.util.concurrent.Executors

import akka.http.scaladsl.Http
import akka.actor.ActorSystem
import akka.http.scaladsl.coding.{Deflate, Gzip, NoCoding}
import akka.http.scaladsl.model.{HttpEntity, _}
import akka.http.scaladsl.unmarshalling.Unmarshal
import akka.stream.{ActorMaterializer, ActorMaterializerSettings}
import centaur.test.metadata.WorkflowMetadata
import centaur.test.workflow.Workflow
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model.headers.HttpEncodings
import akka.util.ByteString
import centaur.api.CromwellStatusJsonSupport._
import centaur.api.CromwellBackendsJsonSupport._
import centaur.{CentaurConfig, SubmittedWorkflow, WorkflowStatus}
import spray.json._

import scala.collection.immutable
import scala.concurrent.{Await, ExecutionContext, Future}
import scala.concurrent.duration.FiniteDuration
import scala.util.Try

object CromwellClient {
  // Do not use scala.concurrent.ExecutionContext.Implicits.global as long as this is using Await.result
  // See https://github.com/akka/akka-http/issues/602
  // And https://github.com/viktorklang/blog/blob/master/Futures-in-Scala-2.12-part-7.md
  final implicit val blockingEc = ExecutionContext.fromExecutor(Executors.newCachedThreadPool())
  
  // Akka HTTP needs both the actor system and a materializer
  final implicit val system = ActorSystem("centaur-acting-like-a-system")
  final implicit val materializer: ActorMaterializer = ActorMaterializer(ActorMaterializerSettings(system))
  final val apiPath = "/api/workflows/v1"

  def decodeResponse(response: HttpResponse): HttpResponse = {
    val decoder = response.encoding match {
      case HttpEncodings.gzip ⇒
        Gzip
      case HttpEncodings.deflate ⇒
        Deflate
      case HttpEncodings.identity ⇒
        NoCoding
    }

    decoder.decode(response)
  }

  def submit(workflow: Workflow): Try[SubmittedWorkflow] = {
    val sourceBodyParts = Map(
      "wdlSource" -> Option(workflow.data.wdl),
      "workflowInputs" -> workflow.data.inputs,
      "workflowOptions" -> insertSecrets(workflow.data.options)
    ) collect { case (name, Some(source)) => Multipart.FormData.BodyPart(name, HttpEntity(MediaTypes.`application/json`, ByteString(source))) }

    val zipBodyParts = Map(
      "wdlDependencies" -> workflow.data.zippedImports
    ) collect { case (name, Some(file)) => Multipart.FormData.BodyPart.fromPath(name, MediaTypes.`application/zip`, file.path) }

    val multipartFormData = Multipart.FormData((sourceBodyParts ++ zipBodyParts).toSeq : _*)
    val requestEntity = multipartFormData.toEntity()

    val submittedWorkflow = for {
      response <- Http().singleRequest(HttpRequest(HttpMethods.POST, CentaurConfig.cromwellUrl + apiPath, immutable.Seq.empty[HttpHeader], requestEntity))
      entity <- response.toEntity.to[CromwellStatus]
    } yield SubmittedWorkflow(UUID.fromString(entity.id), CentaurConfig.cromwellUrl, workflow)
    sendReceiveFutureCompletion(submittedWorkflow)
  }

  def status(workflow: SubmittedWorkflow): Try[WorkflowStatus] = {
    val workflowStatus = for {
      response <- Http().singleRequest(HttpRequest(uri = CentaurConfig.cromwellUrl + apiPath + "/" + workflow.id + "/status"))
      entity <- response.toEntity.to[CromwellStatus]
    } yield WorkflowStatus(entity.status)

    sendReceiveFutureCompletion(workflowStatus)
  }

  def metadata(workflow: SubmittedWorkflow): Try[WorkflowMetadata] = {
    val workflowMetadata = for {
      response <- Http().singleRequest(HttpRequest(uri = CentaurConfig.cromwellUrl + apiPath + "/" + workflow.id + "/metadata"))
      entity <- decodeResponse(response).toEntity.to[String]
    } yield WorkflowMetadata.fromMetadataJson(entity).toOption.get

    sendReceiveFutureCompletion(workflowMetadata)
  }

  lazy val backends: Try[CromwellBackends] = {
    val backends = for {
      response <- Http().singleRequest(HttpRequest(uri = CentaurConfig.cromwellUrl + apiPath + "/backends"))
      entity <- response.toEntity.to[CromwellBackends]
    } yield entity

    sendReceiveFutureCompletion(backends)
  }

  /**
    * Ensure that the Future completes within the specified timeout. If it does not, or if the Future fails,
    * will return a Failure, otherwise a Success
    */
  def awaitFutureCompletion[T](x: Future[T], timeout: FiniteDuration) = Try(Await.result(x, timeout))
  def sendReceiveFutureCompletion[T](x: Future[T]) = awaitFutureCompletion(x, CentaurConfig.sendReceiveTimeout)

  private def insertSecrets(options: Option[String]): Option[String] = {
    import DefaultJsonProtocol._
    val tokenKey = "refresh_token"

    def addToken(optionsMap: Map[String, JsValue]): Map[String, JsValue] = {
      CentaurConfig.optionalToken match {
        case Some(token) if optionsMap.get(tokenKey).isDefined => optionsMap + (tokenKey -> JsString(token))
        case _ => optionsMap
      }
    }

    options map (o => addToken(o.toString.parseJson.asJsObject.convertTo[Map[String, JsValue]]).toJson.toString)
  }

  final implicit class EnhancedHttpResponse(val response: HttpResponse) extends AnyVal {
    /**
      * @throws RuntimeException if the response code is not 2xx, intended to be used inside mapping of Futures
      */
    def toEntity: Unmarshal[ResponseEntity] = response match {
      case HttpResponse(c: StatusCodes.Success, headers, entity, _) => Unmarshal(entity)
      case HttpResponse(code, _, _, _) => throw new RuntimeException(s"$code")
    }
  }
}
