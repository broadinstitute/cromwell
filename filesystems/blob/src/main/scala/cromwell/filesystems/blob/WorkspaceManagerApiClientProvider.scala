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

  def getControlledAzureResourceApi: ControlledAzureResourceApi

}

class HttpWorkspaceManagerClientProvider(baseWorkspaceManagerUrl: WorkspaceManagerURL, token: String) extends WorkspaceManagerApiClientProvider {
  def getApiClient: ApiClient = {
    val client: ApiClient = new ApiClient()
    client.setBasePath(baseWorkspaceManagerUrl.value)
    client.setAccessToken(token)

    client
  }

  def getControlledAzureResourceApi: ControlledAzureResourceApi =
    new ControlledAzureResourceApi(getApiClient)

}