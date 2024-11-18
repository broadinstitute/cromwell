package cromwell.cloudsupport.azure

import com.azure.core.management.AzureEnvironment
import com.typesafe.config.ConfigFactory
import net.ceedubs.ficus.Ficus._

object AzureConfiguration {
  private val conf = ConfigFactory.load().getConfig("azure")
  val azureEnvironment =
    AzureEnvironmentConverter.fromString(
      conf.as[Option[String]]("azure-environment").getOrElse(AzureEnvironmentConverter.Azure)
    )
  val azureTokenScopeManagement = conf.as[String]("token-scope-management")
}

object AzureEnvironmentConverter {
  val Azure: String = "AzureCloud"
  val AzureGov: String = "AzureUSGovernmentCloud"

  def fromString(s: String): AzureEnvironment = s match {
    case AzureGov => AzureEnvironment.AZURE_US_GOVERNMENT
    // a bit redundant, but I want to have a explicit case for Azure for clarity, even though it's the default
    case Azure => AzureEnvironment.AZURE
    case _ => AzureEnvironment.AZURE
  }
}
