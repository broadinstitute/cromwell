package cromwell.services.auth

import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage

object GithubAuthVending {
  sealed trait GithubAuthVendingMessage extends ServiceRegistryMessage {
    override def serviceName: String = "GithubAuthVending"
  }

  // types of tokens
  case class TerraToken(value: String)
  case class GithubToken(value: String)

  case class GithubAuthRequest(terraToken: TerraToken) extends GithubAuthVendingMessage

  sealed trait GithubAuthVendingResponse extends GithubAuthVendingMessage
  case class GithubAuthTokenResponse(githubAccessToken: GithubToken) extends GithubAuthVendingResponse
  case object NoGithubAuthResponse extends GithubAuthVendingResponse
  case class GithubAuthVendingFailure(errorMsg: String) extends GithubAuthVendingResponse

}
