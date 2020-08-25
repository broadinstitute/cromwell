package cromwell.backend.google.pipelines.v2alpha1

import java.net.URL

import cats.effect.IO
import com.google.api.client.http.{HttpRequest, HttpRequestInitializer, HttpResponse}
import com.google.api.services.bigquery.BigqueryScopes
import com.google.api.services.cloudresourcemanager.CloudResourceManager
import com.google.api.services.compute.ComputeScopes
import com.google.api.services.genomics.v2alpha1.model._
import com.google.api.services.genomics.v2alpha1.{Genomics, GenomicsScopes}
import com.google.api.services.oauth2.Oauth2Scopes
import com.google.api.services.storage.StorageScopes
import com.google.auth.http.HttpCredentialsAdapter
import com.typesafe.config.ConfigFactory
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.{GcsTransferConfiguration, VirtualPrivateCloudConfiguration}
import cromwell.backend.google.pipelines.common.action.ActionUtils
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.api.{PipelinesApiFactoryInterface, PipelinesApiRequestFactory}
import cromwell.backend.google.pipelines.common.{GoogleCloudScopes, ProjectLabels}
import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._
import cromwell.backend.google.pipelines.v2alpha1.api._
import cromwell.backend.standard.StandardAsyncJob
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.DockerConfiguration
import cromwell.core.logging.JobLogger
import io.circe.Decoder
import io.circe.generic.semiauto.deriveDecoder
import io.circe.parser.decode
import mouse.all._
import org.apache.commons.lang3.exception.ExceptionUtils
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize

import scala.collection.JavaConverters._

case class GenomicsFactory(applicationName: String, authMode: GoogleAuthMode, endpointUrl: URL)(implicit gcsTransferConfiguration: GcsTransferConfiguration) extends PipelinesApiFactoryInterface
  with ContainerSetup
  with MonitoringAction
  with Localization
  with UserAction
  with Delocalization
  with MemoryRetryCheckAction
  with SSHAccessAction {

  override def build(initializer: HttpRequestInitializer): PipelinesApiRequestFactory = new PipelinesApiRequestFactory {
    implicit lazy val googleProjectMetadataLabelDecoder: Decoder[ProjectLabels] = deriveDecoder

    val ResourceManagerAuthScopes = List(GenomicsScopes.CLOUD_PLATFORM)
    val VirtualPrivateCloudNetworkPath = "projects/%s/global/networks/%s/"

    val genomics: Genomics = new Genomics.Builder(
      GoogleAuthMode.httpTransport,
      GoogleAuthMode.jsonFactory,
      initializer)
      .setApplicationName(applicationName)
      .setRootUrl(endpointUrl.toString)
      .build

    override def cancelRequest(job: StandardAsyncJob): HttpRequest = {
      genomics.projects().operations().cancel(job.jobId, new CancelOperationRequest()).buildHttpRequest()
    }

    override def getRequest(job: StandardAsyncJob): HttpRequest = {
      genomics.projects().operations().get(job.jobId).buildHttpRequest()
    }

    override def runRequest(createPipelineParameters: CreatePipelineParameters, jobLogger: JobLogger): HttpRequest = {

      def projectMetadataRequest(vpcConfig: VirtualPrivateCloudConfiguration): IO[HttpRequest] = {
        IO {
          val workflowOptions = createPipelineParameters.jobDescriptor.workflowDescriptor.workflowOptions
          val credentials = vpcConfig.auth.credentials(workflowOptions.get(_).get, ResourceManagerAuthScopes)

          val httpCredentialsAdapter = new HttpCredentialsAdapter(credentials)
          val cloudResourceManagerBuilder = new CloudResourceManager
          .Builder(GoogleAuthMode.httpTransport, GoogleAuthMode.jsonFactory, httpCredentialsAdapter)
            .setApplicationName(applicationName)
            .build()

          val project = cloudResourceManagerBuilder.projects().get(createPipelineParameters.projectId)

          project.buildHttpRequest()
        }
      }


      def projectMetadataResponseToLabels(httpResponse: HttpResponse): IO[ProjectLabels] = {
        IO.fromEither(decode[ProjectLabels](httpResponse.parseAsString())).handleErrorWith {
          e => IO.raiseError(new RuntimeException(s"Failed to parse labels from project metadata response from Google Cloud Resource Manager API. " +
            s"${ExceptionUtils.getMessage(e)}", e))
        }
      }


      def networkFromLabels(vpcConfig: VirtualPrivateCloudConfiguration, projectLabels: ProjectLabels): Network = {
        val networkLabelOption = projectLabels.labels.find(l => l._1.equals(vpcConfig.name))
        val subnetworkLabelOption = vpcConfig.subnetwork.flatMap(s => projectLabels.labels.find(l => l._1.equals(s)))

        networkLabelOption match {
          case Some(networkLabel) =>
            val network = new Network()
              .setUsePrivateAddress(createPipelineParameters.effectiveNoAddressValue)
              .setName(VirtualPrivateCloudNetworkPath.format(createPipelineParameters.projectId, networkLabel._2))

            subnetworkLabelOption foreach { case(_, subnet) => network.setSubnetwork(subnet) }
            network
          case None =>
            // Falling back to running the job on default network since the project does not provide the custom network
            // specifying keys in its metadata
            new Network().setUsePrivateAddress(createPipelineParameters.effectiveNoAddressValue)
        }
      }


      def createNetworkWithVPC(vpcConfig: VirtualPrivateCloudConfiguration): IO[Network] = {
        for {
          projectMetadataResponse <- projectMetadataRequest(vpcConfig).map(_.executeAsync().get())
          projectLabels <- projectMetadataResponseToLabels(projectMetadataResponse)
          networkLabels <- IO(networkFromLabels(vpcConfig, projectLabels))
        } yield networkLabels
      }


      def createNetwork(): Network = {
        createPipelineParameters.virtualPrivateCloudConfiguration match {
          case None => new Network().setUsePrivateAddress(createPipelineParameters.effectiveNoAddressValue)
          case Some(vpcConfig) =>
            createNetworkWithVPC(vpcConfig).handleErrorWith {
              e => IO.raiseError(new RuntimeException(s"Failed to create Network object for project `${createPipelineParameters.projectId}`. " +
                s"Error(s): ${ExceptionUtils.getMessage(e)}", e))
            }.unsafeRunSync()
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

      val sortedActions =
        ActionUtils.sortActions[Action](
          containerSetup = containerSetup,
          localization = localization,
          userAction = userAction,
          memoryRetryAction = memoryRetryAction,
          deLocalization = deLocalization,
          monitoringSetup = monitoringSetup,
          monitoringShutdown = monitoringShutdown,
          sshAccess = sshAccess,
          isBackground =
            action =>
              Option(action.getFlags)
                .map(_.asScala)
                .toList
                .flatten
                .contains(ActionFlag.RunInBackground.toString),
        )

      val serviceAccount = new ServiceAccount()
        .setEmail(createPipelineParameters.computeServiceAccount)
        .setScopes(
          List(
            GenomicsScopes.GENOMICS,
            ComputeScopes.COMPUTE,
            StorageScopes.DEVSTORAGE_FULL_CONTROL,
            GoogleCloudScopes.KmsScope,
            // Profile and Email scopes are requirements for interacting with Martha v2
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
        .setProjectId(createPipelineParameters.projectId)
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

      genomics.pipelines().run(pipelineRequest).buildHttpRequest()
    }
  }

  override def usesEncryptedDocker: Boolean = true
}

//noinspection ScalaRedundantConversion
object GenomicsFactory {
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
