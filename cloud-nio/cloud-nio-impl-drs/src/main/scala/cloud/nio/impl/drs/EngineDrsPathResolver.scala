package cloud.nio.impl.drs

import com.google.auth.oauth2.{AccessToken, OAuth2Credentials}

import scala.concurrent.duration._


case class EngineDrsPathResolver(drsConfig: DrsConfig,
                                 accessTokenAcceptableTTL: Duration,
                                 authCredentials: OAuth2Credentials,
                                )
  extends DrsPathResolver(drsConfig) {

  //Based on method from GcrRegistry
  override def getAccessToken: String = {
    def accessTokenTTLIsAcceptable(accessToken: AccessToken): Boolean = {
      (accessToken.getExpirationTime.getTime - System.currentTimeMillis()).millis.gteq(accessTokenAcceptableTTL)
    }

    Option(authCredentials.getAccessToken) match {
      case Some(accessToken) if accessTokenTTLIsAcceptable(accessToken) => accessToken.getTokenValue
      case _ =>
        authCredentials.refresh()
        authCredentials.getAccessToken.getTokenValue
    }
  }
}
