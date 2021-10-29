package drs.localizer.tokenproviders

trait AccessTokenProvider {
  def getAccessToken(): String
}
