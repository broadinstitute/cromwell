package cromwell.backend.google.pipelines.v2alpha1

import java.net.URL

import com.google.api.client.http.HttpRequestInitializer
import com.google.api.services.compute.ComputeScopes
import com.google.api.services.genomics.v2alpha1.{Genomics, GenomicsScopes}
import com.google.api.services.genomics.v2alpha1.model._
import com.google.api.services.oauth2.Oauth2Scopes
import com.google.api.services.storage.StorageScopes
import cromwell.backend.google.pipelines.common.PipelinesApiAttributes.LocalizationConfiguration
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.api.{PipelinesApiFactoryInterface, PipelinesApiRequestFactory}
import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._
import cromwell.backend.google.pipelines.v2alpha1.api._
import cromwell.backend.standard.StandardAsyncJob
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.logging.JobLogger
import mouse.all._

import scala.collection.JavaConverters._

case class GenomicsFactory(applicationName: String, authMode: GoogleAuthMode, endpointUrl: URL)(implicit localizationConfiguration: LocalizationConfiguration) extends PipelinesApiFactoryInterface
  with ContainerSetup
  with Localization
  with UserAction
  with Delocalization {

  override def build(initializer: HttpRequestInitializer): PipelinesApiRequestFactory = new PipelinesApiRequestFactory {
    val genomics = new Genomics.Builder(
      GoogleAuthMode.httpTransport,
      GoogleAuthMode.jsonFactory,
      initializer)
      .setApplicationName(applicationName)
      .setRootUrl(endpointUrl.toString)
      .build

    override def cancelRequest(job: StandardAsyncJob) = {
      genomics.projects().operations().cancel(job.jobId, new CancelOperationRequest()).buildHttpRequest()
    }

    override def getRequest(job: StandardAsyncJob) = {
      genomics.projects().operations().get(job.jobId).buildHttpRequest()
    }

    override def runRequest(createPipelineParameters: CreatePipelineParameters, jobLogger: JobLogger) = {
      // Disks defined in the runtime attributes
      val disks = createPipelineParameters |> toDisks
      // Mounts for disks defined in the runtime attributes
      val mounts = createPipelineParameters |> toMounts

      val containerSetup: List[Action] = containerSetupActions(mounts)
      val localization: List[Action] = localizeActions(createPipelineParameters, mounts)
      val userAction: List[Action] = userActions(createPipelineParameters, mounts)
      val deLocalization: List[Action] = deLocalizeActions(createPipelineParameters, mounts)
      val allActions = containerSetup ++ localization ++ userAction ++ deLocalization

      val environment = Map.empty[String, String].asJava

      // Start background actions first, leave the rest as is
      val sortedActions = allActions.sortWith({
        case (a1, _) => Option(a1.getFlags).map(_.asScala).toList.flatten.contains(ActionFlag.RunInBackground.toString)
      })

      val serviceAccount = new ServiceAccount()
        .setEmail(createPipelineParameters.computeServiceAccount)
        .setScopes(
          List(
            GenomicsScopes.GENOMICS,
            ComputeScopes.COMPUTE,
            StorageScopes.DEVSTORAGE_FULL_CONTROL,
            GenomicsFactory.KmsScope,
            // Profile and Email scopes are requirements for interacting with Martha v2
            Oauth2Scopes.USERINFO_EMAIL,
            Oauth2Scopes.USERINFO_PROFILE
          ).asJava
        )

      val network = new Network()
        .setUsePrivateAddress(createPipelineParameters.runtimeAttributes.noAddress)

      val accelerators = createPipelineParameters.runtimeAttributes
        .gpuResource.map(toAccelerator).toList.asJava
      
      /*
       * Adjust for the additional docker images Cromwell uses for (de)localization
       * At the moment, google/cloud-sdk:slim (173MB) and stedolan/jq:latest (182MB)
       * Round it up to 1GB
       */
      val adjustedBootDiskSize = createPipelineParameters.runtimeAttributes.bootDiskSize + 1

      val virtualMachine = new VirtualMachine()
        .setDisks(disks.asJava)
        .setPreemptible(createPipelineParameters.preemptible)
        .setServiceAccount(serviceAccount)
        .setMachineType(createPipelineParameters.runtimeAttributes |> toMachineType(jobLogger))
        .setBootDiskSizeGb(adjustedBootDiskSize)
        .setLabels(createPipelineParameters.labels.asJavaMap)
        .setNetwork(network)
        .setAccelerators(accelerators)

      createPipelineParameters.runtimeAttributes.gpuResource foreach { resource =>
        virtualMachine.setNvidiaDriverVersion(resource.nvidiaDriverVersion)
      }

      val resources = new Resources()
        .setProjectId(createPipelineParameters.projectId)
        .setZones(createPipelineParameters.runtimeAttributes.zones.asJava)
        .setVirtualMachine(virtualMachine)

      val pipeline = new Pipeline()
        .setResources(resources)
        .setActions(sortedActions.asJava)
        .setEnvironment(environment)

      val pipelineRequest = new RunPipelineRequest()
        .setPipeline(pipeline)
        .setLabels(createPipelineParameters.labels.asJavaMap)

      genomics.pipelines().run(pipelineRequest).buildHttpRequest()
    }
  }

  override def usesEncryptedDocker: Boolean = true
}

object GenomicsFactory {
  /**
    * More restricted version of com.google.api.services.cloudkms.v1.CloudKMSScopes.CLOUD_PLATFORM
    * Could use that scope to keep things simple, but docs say to use a more restricted scope:
    *
    *   https://cloud.google.com/kms/docs/accessing-the-api#google_compute_engine
    *
    * For some reason this scope isn't listed as a constant under CloudKMSScopes.
    */
  val KmsScope = "https://www.googleapis.com/auth/cloudkms"
}
