package drs.localizer

import cloud.nio.impl.drs.{DrsConfig, DrsCredentials, DrsPathResolver}
import common.validation.ErrorOr.ErrorOr


class DrsLocalizerDrsPathResolver(drsConfig: DrsConfig, drsCredentials: DrsCredentials) extends DrsPathResolver(drsConfig) {
  override def getAccessToken: ErrorOr[String] = drsCredentials.getAccessToken
}
