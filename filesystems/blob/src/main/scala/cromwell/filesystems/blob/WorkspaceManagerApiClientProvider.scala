package cromwell.filesystems.blob

import bio.terra.workspace.api._
import bio.terra.workspace.client.ApiClient
import bio.terra.workspace.model.{IamRole, ResourceType, StewardshipType}
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
  def getWorkspaceApi(token: String): WsmWorkspaceApi
  def getResourceApi(token: String): WsmResourceApi
  def getLandingZonesApi(token: String): WsmLandingZonesApi
}

class HttpWorkspaceManagerClientProvider(baseWorkspaceManagerUrl: WorkspaceManagerURL) extends WorkspaceManagerApiClientProvider {
  private def getApiClient: ApiClient = {
    val client: ApiClient = new ApiClient()
    client.setBasePath(baseWorkspaceManagerUrl.value)
    client
  }

  def getWorkspaceApi(token: String): WsmWorkspaceApi = {
    val apiClient = getApiClient
    apiClient.setAccessToken(token)
    WsmWorkspaceApi(new WorkspaceApi(apiClient))
  }

  def getResourceApi(token: String): WsmResourceApi = {
    val apiClient = getApiClient
    apiClient.setAccessToken(token)
    WsmResourceApi(new ResourceApi(apiClient))
  }

  def getLandingZonesApi(token: String): WsmLandingZonesApi = {
    val apiClient = getApiClient
    apiClient.setAccessToken(token)
    WsmLandingZonesApi(new LandingZonesApi(apiClient))
  }

  def getControlledAzureResourceApi(token: String): WsmControlledAzureResourceApi = {
    val apiClient = getApiClient
    apiClient.setAccessToken(token)
    WsmControlledAzureResourceApi(new ControlledAzureResourceApi(apiClient))
  }
}

case class WsmLandingZonesApi(landingZonesApi : LandingZonesApi) {
  def findLandingZoneForBillingProfile(billingProfileId: UUID): Try[UUID] = {
    for {
      landingZones <- Try(landingZonesApi.listAzureLandingZones(billingProfileId).getLandingzones())
      landingZoneOption = landingZones.asScala.find(_ => true) // Find the first expecting that there should be only one
      landingZone <- landingZoneOption.map(_.getLandingZoneId()).toRight(new Exception("Billing profile not found")).toTry
    } yield landingZone
  }

  def findStorageAccountForLandingZone(landingZoneId: UUID): Try[StorageAccountName] = {
    for {
      resources <- Try(landingZonesApi.listAzureLandingZoneResources(landingZoneId).getResources())
      sharedResourcesOption = resources.asScala.find(pg => pg.getPurpose() == "SHARED_RESOURCE")
      sharedResources <- sharedResourcesOption.toRight(new Exception("No shared resources found")).toTry
      storageAccountsOption = sharedResources.getDeployedResources().asScala.find(a => a.getResourceType == "Microsoft.Storage/storageAccounts")
      storageAccount <- storageAccountsOption.toRight(new Exception("No storage accounts found")).toTry
      storageAccountId = StorageAccountName(storageAccount.getResourceId())
    } yield storageAccountId
  }
}
case class WsmResourceApi(resourcesApi : ResourceApi) {
  def findContainerResourceId(workspaceId : UUID, container: BlobContainerName): Try[UUID] = {
    for {
      workspaceResources <- Try(resourcesApi.enumerateResources(workspaceId, 0, 1, ResourceType.AZURE_STORAGE_CONTAINER, StewardshipType.CONTROLLED).getResources())
      workspaceStorageContainerOption = workspaceResources.asScala.find(r => r.getMetadata().getName() == container.value)
      workspaceStorageContainer <- workspaceStorageContainerOption.toRight(new Exception("No storage container found for this workspace")).toTry
      resourceId = workspaceStorageContainer.getMetadata().getResourceId()
    } yield resourceId
  }
}
case class WsmWorkspaceApi(workspaceApi : WorkspaceApi) {
  def findBillingProfileId(workspaceId: UUID): Try[UUID] = {
    for {
      wd <- Try(workspaceApi.getWorkspace(workspaceId, IamRole.READER))
      billingProfile = UUID.fromString(wd.getSpendProfile())
    } yield billingProfile
  }
}
case class WsmControlledAzureResourceApi(controlledAzureResourceApi : ControlledAzureResourceApi) {
  def createAzureStorageContainerSasToken(workspaceId: UUID, resourceId: UUID): Try[AzureSasCredential] = {
    for {
      sas <- Try(controlledAzureResourceApi.createAzureStorageContainerSasToken(
        workspaceId,
        resourceId,
        null,
        null,
        null,
        null
      ).getToken)
    } yield new AzureSasCredential(sas)
  }
}
