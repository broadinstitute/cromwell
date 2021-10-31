package drs.localizer.accesstokens

import drs.localizer.CommandLineArguments

case class AzureB2CAccessTokenStrategy(commandLineArguments: CommandLineArguments) extends AccessTokenStrategy {
  // Wire in logic to extract B2C token from KeyVault using UAMI
  override def getAccessToken(): String = ???
}
