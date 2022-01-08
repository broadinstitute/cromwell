package drs.localizer

import cloud.nio.impl.drs.{DrsConfig, DrsPathResolver}
import common.validation.ErrorOr.ErrorOr
import drs.localizer.accesstokens.AccessTokenStrategy


class DrsLocalizerDrsPathResolver(drsConfig: DrsConfig, accessTokenStrategy: AccessTokenStrategy) extends DrsPathResolver(drsConfig) {
  override def getAccessToken: ErrorOr[String] = accessTokenStrategy.getAccessToken()
}
