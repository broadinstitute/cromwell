package drs.localizer.accesstokens

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import drs.localizer.CommandLineArguments

case class AzureB2CAccessTokenStrategy(commandLineArguments: CommandLineArguments) extends AccessTokenStrategy {
  // Wire in logic to extract B2C token from KeyVault using UAMI
  override def getAccessToken(): ErrorOr[String] = "Azure UAMI access token strategy not yet implemented".invalidNel
}
