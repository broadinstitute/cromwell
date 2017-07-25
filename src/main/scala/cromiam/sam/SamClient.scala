package cromiam.sam

import akka.actor.ActorSystem
import akka.http.scaladsl.Http
import akka.http.scaladsl.model.headers.Authorization
import akka.http.scaladsl.model._
import akka.stream.ActorMaterializer
import akka.util.ByteString

import scala.concurrent.{ExecutionContext, ExecutionContextExecutor, Future}
import cats.syntax.functor._
import cats.syntax.traverse._
import cats.instances.future._
import cats.instances.list._
import cromiam.sam.SamClient._
import cromiam.server.config.CromIamServerConfig

trait SamClient {

  implicit def system: ActorSystem
  implicit def executor: ExecutionContextExecutor
  implicit def materializer: ActorMaterializer

  protected def configuration: CromIamServerConfig

  /**
    * Requests authorization for the authenticated user to perform an actions from the Sam service.
    * @return Successful future if the auth is accepted, a Failure otherwise.
    */
  protected def requestAuth(authorizationRequest: WorkflowAuthorizationRequest): Future[Unit] = {
    def validateEntityBytes(byteString: ByteString): Future[Unit] = if (byteString.utf8String == "true") Future.successful(()) else Future.failed(SamDenialException)
    val request = samAuthorizationRequest(authorizationRequest)
    for {
      response <- Http().singleRequest(request)
      entityBytes <- response.entity.dataBytes.runFold(ByteString(""))(_ ++ _)
      _ <- validateEntityBytes(entityBytes)
    } yield ()
  }

  protected def registerCreation(authorization: Authorization, workflowIds: List[String])(implicit ec: ExecutionContext): Future[Unit] = {
    def registerId(id: String) = Http().singleRequest(samRegisterPost(authorization, id)) map { _.discardEntityBytes() }

    (workflowIds traverse registerId).void
  }

  private lazy val SamBaseUri = s"${configuration.samConfig.scheme}://${configuration.samConfig.interface}:${configuration.samConfig.port}/api/resource"
  private def samBaseUriForWorkflow(workflowId: String) = s"$SamBaseUri/workflow/$workflowId"
  private def samRegisterUri(workflowId: String) = Uri(samBaseUriForWorkflow(workflowId))
  private def samAuthorizeActionUri(authorizationRequest: WorkflowAuthorizationRequest) =
    Uri(s"${samBaseUriForWorkflow(authorizationRequest.workflowId)}/action/${authorizationRequest.action}")

  private def samAuthorizationRequest(authorizationRequest: WorkflowAuthorizationRequest): HttpRequest = HttpRequest(
    method = HttpMethods.GET,
    uri = samAuthorizeActionUri(authorizationRequest),
    headers = List[HttpHeader](authorizationRequest.authorization)
  )

  private def samRegisterPost(authorization: Authorization, workflowId: String) = HttpRequest(
    method = HttpMethods.POST,
    uri = samRegisterUri(workflowId),
    headers = List[HttpHeader](authorization)
  )
}

object SamClient {
  case object SamDenialException extends Exception("Access Denied")
  final case class WorkflowAuthorizationRequest(authorization: Authorization, workflowId: String, action: String)
}

