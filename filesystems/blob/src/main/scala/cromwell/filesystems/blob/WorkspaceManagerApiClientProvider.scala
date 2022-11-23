package cromwell.filesystems.blob

import bio.terra.workspace.api.ControlledAzureResourceApi
import bio.terra.workspace.client.ApiClient

/**
  * Represents a way to get various workspace manager clients
  *
  * Pared down from `org.broadinstitute.dsde.rawls.dataaccess.workspacemanager.WorkspaceManagerApiClientProvider`
  *
  * For testing, create an anonymous subclass as in `org.broadinstitute.dsde.rawls.dataaccess.workspacemanager.HttpWorkspaceManagerDAOSpec`
  */
trait WorkspaceManagerApiClientProvider {
  def getApiClient: ApiClient

  def getControlledAzureResourceApi(token: String): ControlledAzureResourceApi
  def getControlledAzureResourceApi(): ControlledAzureResourceApi

}

class HttpWorkspaceManagerClientProvider(baseWorkspaceManagerUrl: WorkspaceManagerURL) extends WorkspaceManagerApiClientProvider {
  def getApiClient: ApiClient = {
    val client: ApiClient = new ApiClient()
    client.setBasePath(baseWorkspaceManagerUrl.value)
    client
  }

  def getControlledAzureResourceApi(token: String): ControlledAzureResourceApi = {
    val apiClient = getApiClient
    apiClient.setAccessToken(token)
    new ControlledAzureResourceApi(apiClient)
  }

  def getControlledAzureResourceApi(): ControlledAzureResourceApi = {
    new ControlledAzureResourceApi(getApiClient)
  }
}
