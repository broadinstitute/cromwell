package cromwell.filesystems.blob

import bio.terra.workspace.api._
import bio.terra.workspace.client.ApiClient
import bio.terra.workspace.model.{ResourceType, StewardshipType}
import com.azure.core.credential.AzureSasCredential

import java.util.UUID
import scala.jdk.CollectionConverters._
import scala.util.Try

/**
  * Represents a way to get a client for interacting with workspace manager controlled resources.
  * Additional WSM clients can be added here if needed.
  *
  * Pared down from `org.broadinstitute.dsde.rawls.dataaccess.workspacemanager.WorkspaceManagerApiClientProvider`
  *
  * For testing, create an anonymous subclass as in `org.broadinstitute.dsde.rawls.dataaccess.workspacemanager.HttpWorkspaceManagerDAOSpec`
  */
trait WorkspaceManagerApiClientProvider {
  def getControlledAzureResourceApi(token: String): WsmControlledAzureResourceApi
  def getResourceApi(token: String): WsmResourceApi
  def getBaseWorkspaceManagerUrl: String
}

class HttpWorkspaceManagerClientProvider(baseWorkspaceManagerUrl: WorkspaceManagerURL)
    extends WorkspaceManagerApiClientProvider {
  private def getApiClient: ApiClient = {
    val client: ApiClient = new ApiClient()
    client.setBasePath(baseWorkspaceManagerUrl.value)
    client
  }

  def getResourceApi(token: String): WsmResourceApi = {
    val apiClient = getApiClient
    apiClient.setAccessToken(token)
    WsmResourceApi(new ResourceApi(apiClient))
  }

  def getControlledAzureResourceApi(token: String): WsmControlledAzureResourceApi = {
    val apiClient = getApiClient
    apiClient.setAccessToken(token)
    WsmControlledAzureResourceApi(new ControlledAzureResourceApi(apiClient))
  }
  def getBaseWorkspaceManagerUrl: String = baseWorkspaceManagerUrl.value
}

case class WsmResourceApi(resourcesApi: ResourceApi) {
  def findContainerResourceId(workspaceId: UUID, container: BlobContainerName): Try[UUID] =
    for {
      workspaceResources <- Try(
        resourcesApi
          .enumerateResources(workspaceId, 0, 10, ResourceType.AZURE_STORAGE_CONTAINER, StewardshipType.CONTROLLED)
          .getResources()
      )
      workspaceStorageContainerOption = workspaceResources.asScala.find(r =>
        r.getMetadata().getName() == container.value
      )
      workspaceStorageContainer <- workspaceStorageContainerOption
        .toRight(new Exception("No storage container found for this workspace"))
        .toTry
      resourceId = workspaceStorageContainer.getMetadata().getResourceId()
    } yield resourceId
}
case class WsmControlledAzureResourceApi(controlledAzureResourceApi: ControlledAzureResourceApi) {
  def createAzureStorageContainerSasToken(workspaceId: UUID, resourceId: UUID): Try[AzureSasCredential] =
    for {
      sas <- Try(
        controlledAzureResourceApi
          .createAzureStorageContainerSasToken(
            workspaceId,
            resourceId,
            null,
            null,
            null,
            null
          )
          .getToken
      )
    } yield new AzureSasCredential(sas)
}
