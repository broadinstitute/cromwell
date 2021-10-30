package drs.localizer

import cloud.nio.impl.drs.{DrsConfig, DrsPathResolver}
import drs.localizer.tokenstrategy.AccessTokenStrategy


class DrsLocalizerDrsPathResolver(drsConfig: DrsConfig, accessTokenProvider: AccessTokenStrategy) extends DrsPathResolver(drsConfig) {
  override def getAccessToken: String = accessTokenProvider.getAccessToken()
}
