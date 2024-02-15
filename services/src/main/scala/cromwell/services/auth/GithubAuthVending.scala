package cromwell.services.auth

import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage

object GithubAuthVending {
  sealed trait GithubAuthVendingMessage extends ServiceRegistryMessage {
    override def serviceName: String = "GithubAuthVending"
  }

  case class GithubAuthRequest(terraToken: String) extends GithubAuthVendingMessage

  sealed trait GithubAuthVendingResponse extends GithubAuthVendingMessage
  case class GithubAuthTokenResponse(accessToken: String) extends GithubAuthVendingResponse
  case object NoGithubAuthResponse extends GithubAuthVendingResponse
  case class GithubAuthVendingFailure(error: Exception) extends GithubAuthVendingResponse

}
