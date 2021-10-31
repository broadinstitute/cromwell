package drs.localizer.accesstokens

trait AccessTokenStrategy {
  def getAccessToken(): String
}
