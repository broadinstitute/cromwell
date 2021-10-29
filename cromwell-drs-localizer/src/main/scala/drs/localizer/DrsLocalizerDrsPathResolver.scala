package drs.localizer

import cloud.nio.impl.drs.{DrsConfig, DrsPathResolver}
import drs.localizer.tokenproviders.AccessTokenProvider


class DrsLocalizerDrsPathResolver(drsConfig: DrsConfig, accessTokenProvider: AccessTokenProvider) extends DrsPathResolver(drsConfig) {
  override def getAccessToken: String = accessTokenProvider.getAccessToken()
}
