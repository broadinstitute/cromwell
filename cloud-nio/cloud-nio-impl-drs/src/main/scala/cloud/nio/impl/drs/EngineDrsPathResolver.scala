package cloud.nio.impl.drs

import common.validation.ErrorOr.ErrorOr

import scala.concurrent.duration._

case class EngineDrsPathResolver(drsConfig: DrsConfig,
                                 accessTokenAcceptableTTL: Duration,
                                 drsCredentials: DrsCredentials,
                                )
  extends DrsPathResolver(drsConfig, retryInternally = false) {

  override def getAccessToken: ErrorOr[String] = {
    drsCredentials.getAccessToken(accessTokenAcceptableTTL)
  }
}
