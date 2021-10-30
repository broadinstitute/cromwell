package drs.localizer.tokenstrategy

trait AccessTokenStrategy {
  def getAccessToken(): String
}
