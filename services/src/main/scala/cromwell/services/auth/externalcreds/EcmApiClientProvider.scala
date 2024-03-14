package cromwell.services.auth.externalcreds

import bio.terra.externalcreds.api.{OauthApi, PublicApi}
import bio.terra.externalcreds.client.ApiClient
import bio.terra.externalcreds.model.{Provider, SystemStatus}
import cats.implicits.catsSyntaxValidatedId
import common.validation.ErrorOr.ErrorOr

import scala.util.{Failure, Success, Try}

trait EcmApiClientProvider {
  def getPublicApi: EcmPublicApi
  def getOauthApi(userToken: String): EcmOauthApi
  def getEcmBaseUrl: String
}

class HttpEcmApiClientProvider(baseEcmUrl: String) extends EcmApiClientProvider {

  override def getPublicApi: EcmPublicApi = {
    val client = new ApiClient()
    client.setBasePath(baseEcmUrl)
    EcmPublicApi(new PublicApi(client))
  }

  override def getOauthApi(userToken: String): EcmOauthApi = {
    val client = new ApiClient()
    client.setBasePath(baseEcmUrl)
    client.setAccessToken(userToken)
    EcmOauthApi(new OauthApi(client))
  }

  override def getEcmBaseUrl: String = baseEcmUrl
}

case class EcmPublicApi(publicApi: PublicApi) {

  def getSystemStatus: SystemStatus =
    publicApi.getStatus
}

case class EcmOauthApi(oauthApi: OauthApi) {
  final private val GithubProvider = Provider.GITHUB

//  def getGithubAccessToken: Try[String] = {
//    Try(
//      oauthApi.getProviderAccessToken(GithubProvider)
//    )
////    match {
////      case Success(token) =>
////        token.validNel
////      case Failure(exception) =>
////        exception.getMessage.invalidNel
////    }
//  }

  def getGithubAccessToken: ErrorOr[String] =
    Try(
      oauthApi.getProviderAccessTokenWithHttpInfo(GithubProvider)
    ) match {
      case Success(token) =>
        token.getBody.validNel
      case Failure(exception) =>
        exception.getMessage.invalidNel
    }
}
