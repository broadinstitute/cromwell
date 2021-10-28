package drs.localizer.tokenproviders

trait AccessTokenProvider {
  def getAccessToken(params: Map[String, String]): String
}
