package cromwell.filesystems.blob

import bio.terra.workspace.api.ControlledAzureResourceApi
import bio.terra.workspace.client.ApiClient

/**
  * Represents a way to get a client for interacting with workspace manager controlled resources.
  * Additional WSM clients can be added here if needed.
  *
  * Pared down from `org.broadinstitute.dsde.rawls.dataaccess.workspacemanager.WorkspaceManagerApiClientProvider`
  *
  * For testing, create an anonymous subclass as in `org.broadinstitute.dsde.rawls.dataaccess.workspacemanager.HttpWorkspaceManagerDAOSpec`
  */
trait WorkspaceManagerApiClientProvider {
  def getControlledAzureResourceApi(token: String): ControlledAzureResourceApi
}

class HttpWorkspaceManagerClientProvider(baseWorkspaceManagerUrl: WorkspaceManagerURL) extends WorkspaceManagerApiClientProvider {
  private def getApiClient: ApiClient = {
    val client: ApiClient = new ApiClient()
    client.setBasePath(baseWorkspaceManagerUrl.value)
    client
  }

  def getControlledAzureResourceApi(token: String): ControlledAzureResourceApi = {
    val apiClient = getApiClient
    apiClient.setAccessToken(token)
    new ControlledAzureResourceApi(apiClient)
  }
}
