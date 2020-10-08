package cromwell.backend.google.pipelines.v2beta

import java.net.URL

import com.google.api.client.http.{HttpRequest, HttpRequestInitializer}
import com.google.api.services.bigquery.BigqueryScopes
import com.google.api.services.compute.ComputeScopes
import com.google.api.services.lifesciences.v2beta.model._
import com.google.api.services.lifesciences.v2beta.CloudLifeSciences
import com.google.api.services.oauth2.Oauth2Scopes
import com.google.api.services.storage.StorageScopes
import com.typesafe.config.ConfigFactory
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.pipelines.common.action.ActionUtils
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.api.{PipelinesApiFactoryInterface, PipelinesApiRequestFactory}
import cromwell.backend.google.pipelines.common.{GoogleCloudScopes, VpcAndSubnetworkProjectLabelValues}
import cromwell.backend.google.pipelines.v2beta.PipelinesConversions._
import cromwell.backend.google.pipelines.v2beta.api._
import cromwell.backend.standard.StandardAsyncJob
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.DockerConfiguration
import cromwell.core.logging.JobLogger
import mouse.all._
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize

import scala.collection.JavaConverters._

case class LifeSciencesFactory(applicationName: String, authMode: GoogleAuthMode, endpointUrl: URL, location: String)(implicit gcsTransferConfiguration: GcsTransferConfiguration) extends PipelinesApiFactoryInterface
  with ContainerSetup
  with MonitoringAction
  with Localization
  with UserAction
  with Delocalization
  with MemoryRetryCheckAction
  with SSHAccessAction {

  override def build(initializer: HttpRequestInitializer): PipelinesApiRequestFactory = new PipelinesApiRequestFactory {
    val VirtualPrivateCloudNetworkPath = "projects/%s/global/networks/%s/"

    val lifeSciences: CloudLifeSciences = new CloudLifeSciences.Builder(
      GoogleAuthMode.httpTransport,
      GoogleAuthMode.jsonFactory,
      initializer)
      .setApplicationName(applicationName)
      .setRootUrl(endpointUrl.toString)
      .build

    override def cancelRequest(job: StandardAsyncJob): HttpRequest = {
      lifeSciences.projects().locations().operations().cancel(job.jobId, new CancelOperationRequest()).buildHttpRequest()
    }

    override def getRequest(job: StandardAsyncJob): HttpRequest = {
      lifeSciences.projects().locations().operations().get(job.jobId).buildHttpRequest()
    }

    override def runRequest(createPipelineParameters: CreatePipelineParameters, jobLogger: JobLogger): HttpRequest = {
      def createNetworkWithVPC(vpcAndSubnetworkProjectLabelValues: VpcAndSubnetworkProjectLabelValues): Network = {
        val network = new Network()
          .setUsePrivateAddress(createPipelineParameters.effectiveNoAddressValue)
          .setNetwork(VirtualPrivateCloudNetworkPath.format(createPipelineParameters.projectId, vpcAndSubnetworkProjectLabelValues.vpcName))

        vpcAndSubnetworkProjectLabelValues.subnetNameOpt.foreach(subnet => network.setSubnetwork(subnet))
        network
      }

      def createNetwork(): Network = {
        createPipelineParameters.vpcNetworkAndSubnetworkProjectLabels match {
          case Some(vpcAndSubnetworkProjectLabelValues) => createNetworkWithVPC(vpcAndSubnetworkProjectLabelValues)
          case _ => new Network().setUsePrivateAddress(createPipelineParameters.effectiveNoAddressValue)
        }
      }

      val allDisksToBeMounted = createPipelineParameters.adjustedSizeDisks ++ createPipelineParameters.referenceDisksForLocalization

      // Disks defined in the runtime attributes and reference-files-localization disks
      val disks = allDisksToBeMounted |> toDisks
      // Mounts for disks defined in the runtime attributes and reference-files-localization disks
      val mounts = allDisksToBeMounted |> toMounts

      val containerSetup: List[Action] = containerSetupActions(mounts)
      val localization: List[Action] = localizeActions(createPipelineParameters, mounts)
      val userAction: List[Action] = userActions(createPipelineParameters, mounts)
      val memoryRetryAction: List[Action] = checkForMemoryRetryActions(createPipelineParameters, mounts)
      val deLocalization: List[Action] = deLocalizeActions(createPipelineParameters, mounts)
      val monitoringSetup: List[Action] = monitoringSetupActions(createPipelineParameters, mounts)
      val monitoringShutdown: List[Action] = monitoringShutdownActions(createPipelineParameters)
      val sshAccess: List[Action] = sshAccessActions(createPipelineParameters, mounts)

      // adding memory as environment variables makes it easy for a user to retrieve the new value of memory
      // on the machine to utilize in their command blocks if needed
      val runtimeMemory = createPipelineParameters.runtimeAttributes.memory
      val environment = Map("MEM_UNIT" -> runtimeMemory.unit.toString, "MEM_SIZE" -> runtimeMemory.amount.toString).asJava

      val sortedActions: List[Action] =
        ActionUtils.sortActions[Action](
          containerSetup = containerSetup,
          localization = localization,
          userAction = userAction,
          memoryRetryAction = memoryRetryAction,
          deLocalization = deLocalization,
          monitoringSetup = monitoringSetup,
          monitoringShutdown = monitoringShutdown,
          sshAccess = sshAccess,
          isBackground = _.getRunInBackground,
        )

      val serviceAccount = new ServiceAccount()
        .setEmail(createPipelineParameters.computeServiceAccount)
        .setScopes(
          List(
            ComputeScopes.COMPUTE,
            StorageScopes.DEVSTORAGE_FULL_CONTROL,
            GoogleCloudScopes.KmsScope,
            // Profile and Email scopes are requirements for interacting with Martha
            Oauth2Scopes.USERINFO_EMAIL,
            Oauth2Scopes.USERINFO_PROFILE,
            // Monitoring scope as POC
            GoogleCloudScopes.MonitoringWrite,
            // Allow read/write with BigQuery
            BigqueryScopes.BIGQUERY
          ).asJava
        )

      val network: Network = createNetwork()

      val accelerators = createPipelineParameters.runtimeAttributes
        .gpuResource.map(toAccelerator).toList.asJava

      /*
       * Adjust using docker images used by Cromwell as well as the tool's docker image size if available
       */
      val adjustedBootDiskSize = {
        val fromRuntimeAttributes = createPipelineParameters.runtimeAttributes.bootDiskSize
        // Compute the decompressed size based on the information available
        val userCommandImageSizeInBytes = createPipelineParameters.jobDescriptor.dockerSize.map(_.toFullSize(DockerConfiguration.instance.sizeCompressionFactor)).getOrElse(0L)
        val userCommandImageSizeInGB = MemorySize(userCommandImageSizeInBytes.toDouble, MemoryUnit.Bytes).to(MemoryUnit.GB).amount
        val userCommandImageSizeRoundedUpInGB = userCommandImageSizeInGB.ceil.toInt

        val totalSize = fromRuntimeAttributes + userCommandImageSizeRoundedUpInGB + ActionUtils.cromwellImagesSizeRoundedUpInGB
        jobLogger.info(s"Adjusting boot disk size to $totalSize GB: $fromRuntimeAttributes GB (runtime attributes) + $userCommandImageSizeRoundedUpInGB GB (user command image) + ${ActionUtils.cromwellImagesSizeRoundedUpInGB} GB (Cromwell support images)")

        totalSize
      }

      val virtualMachine = new VirtualMachine()
        .setDisks(disks.asJava)
        .setPreemptible(createPipelineParameters.preemptible)
        .setServiceAccount(serviceAccount)
        .setMachineType(createPipelineParameters.runtimeAttributes |> toMachineType(jobLogger))
        .setBootDiskSizeGb(adjustedBootDiskSize)
        .setLabels(createPipelineParameters.googleLabels.map(label => label.key -> label.value).toMap.asJava)
        .setNetwork(network)
        .setAccelerators(accelerators)

      createPipelineParameters.runtimeAttributes.gpuResource foreach { resource =>
        virtualMachine.setNvidiaDriverVersion(resource.nvidiaDriverVersion)
      }

      createPipelineParameters.runtimeAttributes.cpuPlatform.map(virtualMachine.setCpuPlatform)

      val resources = new Resources()
        .setZones(createPipelineParameters.runtimeAttributes.zones.asJava)
        .setVirtualMachine(virtualMachine)

      val pipeline = new Pipeline()
        .setResources(resources)
        .setActions(sortedActions.asJava)
        .setEnvironment(environment)
        .setTimeout(createPipelineParameters.pipelineTimeout.toSeconds + "s")

      val pipelineRequest = new RunPipelineRequest()
        .setPipeline(pipeline)
        .setLabels(createPipelineParameters.googleLabels.map(label => label.key -> label.value).toMap.asJava)

      val parent = s"projects/${createPipelineParameters.projectId}/locations/$location"
      lifeSciences.projects().locations().pipelines().run(parent, pipelineRequest).buildHttpRequest()
    }
  }

  override def usesEncryptedDocker: Boolean = true
}

//noinspection ScalaRedundantConversion
object LifeSciencesFactory {
  private val config = ConfigFactory.load().getConfig("google")

  /**
    * An image with the Google Cloud SDK installed.
    * http://gcr.io/google.com/cloudsdktool/cloud-sdk
    *
    * FYI additional older versions are available on DockerHub at:
    * https://hub.docker.com/r/google/cloud-sdk
    *
    * When updating this value, also consider updating the CromwellImagesSizeRoundedUpInGB below.
    */
  val CloudSdkImage: String = if (config.hasPath("cloud-sdk-image-url")) { config.getString("cloud-sdk-image-url").toString } else "gcr.io/google.com/cloudsdktool/cloud-sdk:276.0.0-slim"

  /*
   * At the moment, cloud-sdk (924MB for 276.0.0-slim) and stedolan/jq (182MB) decompressed ~= 1.1 GB
   */
  val CromwellImagesSizeRoundedUpInGB: Int = if (config.hasPath("cloud-sdk-image-size-gb")) { config.getInt("cloud-sdk-image-size-gb") } else 1
}
