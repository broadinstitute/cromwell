package drs.localizer.accesstokens

import cats.syntax.validated._
import com.azure.core.credential.TokenRequestContext
import com.azure.identity.DefaultAzureCredentialBuilder
import common.validation.ErrorOr.ErrorOr
import drs.localizer.CommandLineArguments

case class AzureAccessTokenStrategy(commandLineArguments: CommandLineArguments) extends AccessTokenStrategy {
  override def getAccessToken(): ErrorOr[String] = {
    val credentials = commandLineArguments.azureIdentityClientId.foldLeft(new DefaultAzureCredentialBuilder()) {
      (builder, clientId) => builder.managedIdentityClientId(clientId)
    }.build()
    val tokenRequestContext = new TokenRequestContext()
    tokenRequestContext.addScopes(".default")
    credentials.getToken(tokenRequestContext).block().getToken.validNel
  }
}
