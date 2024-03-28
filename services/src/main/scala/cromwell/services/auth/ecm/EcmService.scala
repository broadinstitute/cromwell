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
     ECM does generally return standard JSON error response, but for 401 status code it seems some other layer in
     between (like the apache proxies, etc) returns HTML pages. This helper method returns custom error message for 401
     status code as it contains HTML tags. For all other status code, the response format is generally of ErrorReport
     schema and this method tries to extract the actual message from the JSON object and return it. In case it fails
     to parse JSON, it returns the original error response body.
     ErrorReport schema: {"message":"<actual_error_msg>", "statusCode":<code>}
   */
  def extractErrorMessage(errorCode: StatusCode, responseBodyAsStr: String): String =
    errorCode match {
      case StatusCodes.Unauthorized => "Invalid or missing authentication credentials."
      case _ =>
        Try(responseBodyAsStr.parseJson) match {
          case Success(JsObject(fields)) =>
            fields.get("message").map(_.toString().replaceAll("\"", "")).getOrElse(responseBodyAsStr)
          case _ => responseBodyAsStr
        }
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
