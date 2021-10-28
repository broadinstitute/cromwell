package drs.localizer.tokenproviders

class AzureB2CTokenProvider extends AccessTokenProvider {
  // Wire up logic to extract B2C token from KeyVault using UAMI
  override def getAccessToken(params: Map[String, String]): String = ???
}
