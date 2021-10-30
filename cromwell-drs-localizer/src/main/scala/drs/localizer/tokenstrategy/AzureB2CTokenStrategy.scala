package drs.localizer.tokenstrategy

import drs.localizer.CommandLineArguments

case class AzureB2CTokenStrategy(commandLineArguments: CommandLineArguments) extends AccessTokenStrategy {
  // Wire up logic to extract B2C token from KeyVault using UAMI
  override def getAccessToken(): String = ???
}
