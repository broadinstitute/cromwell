package cromwell.services.auth.ecm

import akka.actor.ActorSystem
import akka.http.scaladsl.Http
import akka.http.scaladsl.model._
import akka.http.scaladsl.model.headers.RawHeader
import akka.util.ByteString

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Success, Try}
import spray.json._

class EcmService(baseEcmUrl: String) {
  private val getGithubAccessTokenApiPath = "api/oauth/v1/github/access-token"

  // TODO: add comments
  private def extractErrorMessage(errorCode: StatusCode, responseBodyAsStr: String): String =
    errorCode match {
      case StatusCodes.Unauthorized =>
        "401 Unauthorized. Invalid or missing authentication credentials."
      case StatusCodes.Forbidden =>
        "403 Forbidden. User doesn't have the right permission(s) to fetch Github token."
      case StatusCodes.BadRequest | StatusCodes.NotFound | StatusCodes.InternalServerError =>
        Try(responseBodyAsStr.parseJson) match {
          case Success(JsObject(fields)) =>
            fields
              .get("message")
              .map(v => s"${errorCode.value}. ${v.toString()}")
              .getOrElse(responseBodyAsStr)
          case _ =>
            responseBodyAsStr
        }
      case _ => responseBodyAsStr
    }

  def getGithubAccessToken(
    userToken: String
  )(implicit actorSystem: ActorSystem, ec: ExecutionContext): Future[String] = {

    def responseEntityToFutureStr(responseEntity: ResponseEntity): Future[String] =
      responseEntity.dataBytes.runFold(ByteString(""))(_ ++ _).map(_.utf8String)

    val headers: HttpHeader = RawHeader("Authorization", s"Bearer $userToken")
    val httpRequest =
      HttpRequest(method = HttpMethods.GET, uri = s"$baseEcmUrl/$getGithubAccessTokenApiPath").withHeaders(headers)

    Http()
      .singleRequest(httpRequest)
      .flatMap((response: HttpResponse) =>
        if (response.status.isFailure()) {
          responseEntityToFutureStr(response.entity) flatMap { errorBody =>
            val errorMessage = extractErrorMessage(response.status, errorBody)
            Future.failed(
              new RuntimeException(s"HTTP $errorMessage")
            )
          }
        } else responseEntityToFutureStr(response.entity)
      )
  }
}
