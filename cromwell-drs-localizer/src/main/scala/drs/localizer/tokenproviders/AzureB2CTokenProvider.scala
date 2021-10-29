package drs.localizer.tokenproviders

import drs.localizer.CommandLineArguments

case class AzureB2CTokenProvider(commandLineArguments: CommandLineArguments) extends AccessTokenProvider {
  // Wire up logic to extract B2C token from KeyVault using UAMI
  override def getAccessToken(): String = ???
}
