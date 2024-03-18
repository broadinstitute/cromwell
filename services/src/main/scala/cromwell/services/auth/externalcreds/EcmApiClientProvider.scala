package cromwell.services.auth.externalcreds

import bio.terra.externalcreds.api.OauthApi
import bio.terra.externalcreds.client.ApiClient
import bio.terra.externalcreds.model.Provider
import cats.implicits.catsSyntaxValidatedId
import common.validation.ErrorOr.ErrorOr
import org.springframework.http.HttpStatus
import org.springframework.web.client.HttpClientErrorException

import scala.util.{Failure, Success, Try}
import spray.json._

trait EcmApiClientProvider {
  def getOauthApi(userToken: String): EcmOauthApi
  def getEcmBaseUrl: String
}

class HttpEcmApiClientProvider(baseEcmUrl: String) extends EcmApiClientProvider {

  override def getOauthApi(userToken: String): EcmOauthApi = {
    val client = new ApiClient()
    client.setBasePath(baseEcmUrl)
    client.setAccessToken(userToken)
    EcmOauthApi(new OauthApi(client))
  }

  override def getEcmBaseUrl: String = baseEcmUrl
}

case class EcmOauthApi(oauthApi: OauthApi) {

  /*
   ECM client returns a RestClientException and error messages for different response codes don't have a standard format.
   This helper method tries to extract the actual error message embedded inside it and returns it. If it fails to
   extract it the message contained in the exception is returned.
   */
  private def extractErrorMessage(exception: Throwable): String = {

    // if the exception message format in the case it ErrorReport JSON in it, is generally as below:
    //    "<status_code> <status_code_test>: {<EOL> "message":<actual_error_msg>,<EOL> "statusCode":<code><EOL>}"
    // this method attempts to extract the value associated with "message" key and return it. In case its not
    // successful it will return exception.getMessage()
    def getMessageFromErrorReport(e: HttpClientErrorException) = {
      Try {
        // remove any HTML tags or <EOL> tags from string
        val strWithoutTags = exception.getMessage.replaceAll("<[^>]*>", "")
        // find the JSON object in the string
        val strBetweenBracesOpt = """\{(.*?)}""".r.findFirstIn(strWithoutTags)
        // if JSON object is found in the string, parse it as JSON and return the error message using "message" key
        // otherwise return the message found in exception
        strBetweenBracesOpt match {
          case Some(jsonAsStr) =>
            Try(jsonAsStr.parseJson) match {
              case Success(JsObject(fields)) =>
                fields.get("message").map(v => s"${e.getStatusCode.value()} ${e.getStatusText}. ${v.toString().replaceAll("\"", "")}").getOrElse(exception.getMessage)
              case _ =>
                exception.getMessage
            }
          case None =>
            exception.getMessage
        }
      } match {
        case Success(errorMsg) => errorMsg
        case Failure(_) => exception.getMessage
      }
    }

    // this returns a generic error message in case the response is 401 or 403 as the exception message returned from
    // ECM in these cases contains HTML tags. In case of 400, 404 or 500 response ECM client seems to be returning error
    // messages in form of ErrorReport schema, in these cases it will try to extract the embedded error message
    // and return it
    exception match {
      case e: HttpClientErrorException => e.getStatusCode match {
        case HttpStatus.UNAUTHORIZED =>
          "401 Unauthorized. Invalid or missing authentication credentials."
        case HttpStatus.FORBIDDEN =>
          "403 Forbidden. User doesn't have the right permission(s) to fetch Github token."
        case HttpStatus.BAD_REQUEST | HttpStatus.NOT_FOUND | HttpStatus.INTERNAL_SERVER_ERROR =>
          getMessageFromErrorReport(e)
        case _ =>
          exception.getMessage
      }
      case _ =>
        exception.getMessage
    }
  }

  def getGithubAccessToken: ErrorOr[String] =
    Try(
      oauthApi.getProviderAccessToken(Provider.GITHUB)
    ) match {
      case Success(token) =>
        token.validNel
      case Failure(exception) =>
        extractErrorMessage(exception).invalidNel
    }
}
