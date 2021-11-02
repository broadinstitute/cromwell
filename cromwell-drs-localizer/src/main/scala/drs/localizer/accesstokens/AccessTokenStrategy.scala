package drs.localizer.accesstokens

import common.validation.ErrorOr.ErrorOr

trait AccessTokenStrategy {
  def getAccessToken(): ErrorOr[String]
}
