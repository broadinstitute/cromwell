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

  /*
     ECM doesn't have a standard error response format. Some of the responses contains HTML tags in it. This helper
     method returns custom error message for 401 and 403 errors as they contain HTML tags. For 400, 404 and 500 the
     Swagger suggests that the response format is of ErrorReport schema and this method tries to extract the
     actual message from the JSON object and returns it. In case of other status codes or if it fails to parse JSON it
     returns the original error response.
     ErrorReport schema: {"message":<actual_error_msg>, "statusCode":<code>}
   */
  def extractErrorMessage(errorCode: StatusCode, responseBodyAsStr: String): String =
    errorCode match {
      case StatusCodes.Unauthorized => "Invalid or missing authentication credentials."
      case StatusCodes.Forbidden => "User doesn't have the right permission(s) to fetch Github token."
      case StatusCodes.BadRequest | StatusCodes.NotFound | StatusCodes.InternalServerError =>
        Try(responseBodyAsStr.parseJson) match {
          case Success(JsObject(fields)) =>
            fields.get("message").map(_.toString().replaceAll("\"", "")).getOrElse(responseBodyAsStr)
          case _ => responseBodyAsStr
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
            Future.failed(new RuntimeException(s"HTTP ${response.status.value}. $errorMessage"))
          }
        } else responseEntityToFutureStr(response.entity)
      )
  }
}
