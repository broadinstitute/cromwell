package cloud.nio.impl.drs

import com.google.auth.oauth2.{AccessToken, OAuth2Credentials}
import common.validation.ErrorOr.ErrorOr

import scala.concurrent.duration._
import cats.syntax.validated._


case class EngineDrsPathResolver(drsConfig: DrsConfig,
                                 accessTokenAcceptableTTL: Duration,
                                 authCredentials: OAuth2Credentials,
                                )
  extends DrsPathResolver(drsConfig, retryInternally = false) {

  //Based on method from GoogleRegistry
  override def getAccessToken: ErrorOr[String] = {
    def accessTokenTTLIsAcceptable(accessToken: AccessToken): Boolean = {
      (accessToken.getExpirationTime.getTime - System.currentTimeMillis()).millis.gteq(accessTokenAcceptableTTL)
    }

    Option(authCredentials.getAccessToken) match {
      case Some(accessToken) if accessTokenTTLIsAcceptable(accessToken) =>
        accessToken.getTokenValue.validNel
      case _ =>
        authCredentials.refresh()
        Option(authCredentials.getAccessToken.getTokenValue) match {
          case Some(accessToken) => accessToken.validNel
          case None => "Could not refresh access token".invalidNel
        }
    }
  }
}
