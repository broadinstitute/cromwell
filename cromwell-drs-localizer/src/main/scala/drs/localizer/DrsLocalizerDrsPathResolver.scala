package drs.localizer

import cloud.nio.impl.drs.{DrsConfig, DrsPathResolver}
import drs.localizer.accesstokens.AccessTokenStrategy


class DrsLocalizerDrsPathResolver(drsConfig: DrsConfig, accessTokenStrategy: AccessTokenStrategy) extends DrsPathResolver(drsConfig) {
  override def getAccessToken: String = accessTokenStrategy.getAccessToken()
}
