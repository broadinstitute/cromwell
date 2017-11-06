package cromiam.sam

import akka.actor.ActorSystem
import akka.http.scaladsl.Http
import akka.http.scaladsl.model.headers.Authorization
import akka.http.scaladsl.model._
import akka.stream.ActorMaterializer
import akka.util.ByteString
import scala.concurrent.{ExecutionContext, ExecutionContextExecutor, Future}
import cats.syntax.traverse._
import cats.instances.list._
import cats.syntax.functor._
import cats.instances.future._
import com.softwaremill.sttp._
import cromiam.sam.SamClient._
import cromiam.server.status.StatusCheckedSubsystem

class SamClient(scheme: String, interface: String, port: Int)(implicit system: ActorSystem,
                                                              ece: ExecutionContextExecutor,
                                                              materializer: ActorMaterializer) extends StatusCheckedSubsystem {
  override val statusUri = uri"$samBaseUri/status"

  /**
    * Requests authorization for the authenticated user to perform an actions from the Sam service.
    * @return Successful future if the auth is accepted, a Failure otherwise.
    */
  def requestAuth(authorizationRequest: WorkflowAuthorizationRequest): Future[Unit] = {
    def validateEntityBytes(byteString: ByteString): Future[Unit] = if (byteString.utf8String == "true") Future.successful(()) else Future.failed(SamDenialException)
    val request = samAuthorizationRequest(authorizationRequest)
    for {
      response <- Http().singleRequest(request)
      entityBytes <- response.entity.dataBytes.runFold(ByteString(""))(_ ++ _)
      _ <- validateEntityBytes(entityBytes)
    } yield ()
  }

  def registerCreation(authorization: Authorization, workflowIds: List[String])(implicit ec: ExecutionContext): Future[Unit] = {
    def registerId(id: String) = Http().singleRequest(samRegisterPost(authorization, id)) map { _.discardEntityBytes() }

    (workflowIds traverse registerId).void
  }

  private lazy val samBaseUri = s"$scheme://$interface:$port"
  private lazy val samBaseResourceUri = s"$samBaseUri/api/resource"

  private def samBaseUriForWorkflow(workflowId: String) = s"$samBaseResourceUri/workflow-collection/$workflowId"
  private def samRegisterUri(workflowId: String) = akka.http.scaladsl.model.Uri(samBaseUriForWorkflow(workflowId))
  private def samAuthorizeActionUri(authorizationRequest: WorkflowAuthorizationRequest) = {
    akka.http.scaladsl.model.Uri(s"${samBaseUriForWorkflow(authorizationRequest.workflowId)}/action/${authorizationRequest.action}")
  }

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

