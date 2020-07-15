package cromwell.backend.google.pipelines.v2beta

import java.net.URL

import cats.effect.IO
import com.google.api.client.http.{HttpRequest, HttpRequestInitializer, HttpResponse}
import com.google.api.services.bigquery.BigqueryScopes
import com.google.api.services.cloudresourcemanager.CloudResourceManager
import com.google.api.services.compute.ComputeScopes
import com.google.api.services.lifesciences.v2beta.model._
import com.google.api.services.lifesciences.v2beta.{CloudLifeSciences, CloudLifeSciencesScopes}
import com.google.api.services.oauth2.Oauth2Scopes
import com.google.api.services.storage.StorageScopes
import com.google.auth.http.HttpCredentialsAdapter
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.{GcsTransferConfiguration, VirtualPrivateCloudConfiguration}
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.api.{PipelinesApiFactoryInterface, PipelinesApiRequestFactory}
import cromwell.backend.google.pipelines.v2beta.PipelinesConversions._
import cromwell.backend.google.pipelines.v2beta.api._
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

case class LifeSciencesFactory(applicationName: String, authMode: GoogleAuthMode, endpointUrl: URL, location: String)(implicit gcsTransferConfiguration: GcsTransferConfiguration) extends PipelinesApiFactoryInterface
  with ContainerSetup
  with MonitoringAction
  with Localization
  with UserAction
  with Delocalization
  with MemoryRetryCheckAction
  with SSHAccessAction {

  override def build(initializer: HttpRequestInitializer): PipelinesApiRequestFactory = new PipelinesApiRequestFactory {
    implicit lazy val googleProjectMetadataLabelDecoder: Decoder[ProjectLabels] = deriveDecoder

    val ResourceManagerAuthScopes = List(CloudLifeSciencesScopes.CLOUD_PLATFORM)
    val VirtualPrivateCloudNetworkPath = "projects/%s/global/networks/%s/"

    val lifeSciences = new CloudLifeSciences.Builder(
      GoogleAuthMode.httpTransport,
      GoogleAuthMode.jsonFactory,
      initializer)
      .setApplicationName(applicationName)
      .setRootUrl(endpointUrl.toString)
      .build

    override def cancelRequest(job: StandardAsyncJob) = {
      lifeSciences.projects().locations().operations().cancel(job.jobId, new CancelOperationRequest()).buildHttpRequest()
    }

    override def getRequest(job: StandardAsyncJob) = {
      lifeSciences.projects().locations().operations().get(job.jobId).buildHttpRequest()
    }

    override def runRequest(createPipelineParameters: CreatePipelineParameters, jobLogger: JobLogger): HttpRequest = {

      def projectMetadataRequest(vpcConfig: VirtualPrivateCloudConfiguration): IO[HttpRequest] = {
        IO {
          val workflowOptions = createPipelineParameters.jobDescriptor.workflowDescriptor.workflowOptions
          val credentials = vpcConfig.auth.credentials(workflowOptions.get(_).get, ResourceManagerAuthScopes)

          val httpCredentialsAdapter = new HttpCredentialsAdapter(credentials)
          val cloudResourceManagerBuilder = new CloudResourceManager
          .Builder(GoogleAuthMode.httpTransport, GoogleAuthMode.jsonFactory, httpCredentialsAdapter)
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
              .setNetwork(VirtualPrivateCloudNetworkPath.format(createPipelineParameters.projectId, networkLabel._2))

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

      // Disks defined in the runtime attributes
      val disks = createPipelineParameters |> toDisks
      // Mounts for disks defined in the runtime attributes
      val mounts = createPipelineParameters |> toMounts

      val containerSetup: List[Action] = containerSetupActions(mounts)
      val localization: List[Action] = localizeActions(createPipelineParameters, mounts)
      val userAction: List[Action] = userActions(createPipelineParameters, mounts)
      val memoryRetryAction: List[Action] = checkForMemoryRetryActions(createPipelineParameters, mounts)
      val deLocalization: List[Action] = deLocalizeActions(createPipelineParameters, mounts)
      val monitoring: List[Action] = monitoringActions(createPipelineParameters, mounts)
      val sshAccess: List[Action] = sshAccessActions(createPipelineParameters, mounts)
      val allActions = containerSetup ++ localization ++ userAction ++ memoryRetryAction ++ deLocalization ++ monitoring ++ sshAccess

      // adding memory as environment variables makes it easy for a user to retrieve the new value of memory
      // on the machine to utilize in their command blocks if needed
      val runtimeMemory = createPipelineParameters.runtimeAttributes.memory
      val environment = Map("MEM_UNIT" -> runtimeMemory.unit.toString, "MEM_SIZE" -> runtimeMemory.amount.toString).asJava

      // Start background actions first, leave the rest as is
      val sortedActions = allActions.sortWith({
        case (a1, _) => a1.getRunInBackground
      })

      val serviceAccount = new ServiceAccount()
        .setEmail(createPipelineParameters.computeServiceAccount)
        .setScopes(
          List(
            ComputeScopes.COMPUTE,
            StorageScopes.DEVSTORAGE_FULL_CONTROL,
            LifeSciencesFactory.KmsScope,
            // Profile and Email scopes are requirements for interacting with Martha v2
            Oauth2Scopes.USERINFO_EMAIL,
            Oauth2Scopes.USERINFO_PROFILE,
            // Monitoring scope as POC
            LifeSciencesFactory.MonitoringWrite,
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
        val actualSizeInBytes = createPipelineParameters.jobDescriptor.dockerSize.map(_.toFullSize(DockerConfiguration.instance.sizeCompressionFactor)).getOrElse(0L)
        val actualSizeInGB = MemorySize(actualSizeInBytes.toDouble, MemoryUnit.Bytes).to(MemoryUnit.GB).amount
        val actualSizeRoundedUpInGB = actualSizeInGB.ceil.toInt

        // If the size of the image is larger than what is in the runtime attributes, adjust it to uncompressed size
        if (actualSizeRoundedUpInGB > fromRuntimeAttributes) {
          jobLogger.info(s"Adjusting boot disk size to $actualSizeRoundedUpInGB GB to account for docker image size")
        }

        Math.max(fromRuntimeAttributes, actualSizeRoundedUpInGB) + LifeSciencesFactory.CromwellImagesSizeRoundedUpInGB
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
        .setTimeout(createPipelineParameters.pipelineTimeout.toSeconds.toString + "s")

      val pipelineRequest = new RunPipelineRequest()
        .setPipeline(pipeline)
        .setLabels(createPipelineParameters.googleLabels.map(label => label.key -> label.value).toMap.asJava)

      val parent = s"projects/${createPipelineParameters.projectId}/locations/$location"
      lifeSciences.projects().locations().pipelines().run(parent, pipelineRequest).buildHttpRequest()
    }
  }

  override def usesEncryptedDocker: Boolean = true
}

object LifeSciencesFactory {
  /**
    * More restricted version of com.google.api.services.cloudkms.v1.CloudKMSScopes.CLOUD_PLATFORM
    * Could use that scope to keep things simple, but docs say to use a more restricted scope:
    *
    *   https://cloud.google.com/kms/docs/accessing-the-api#google_compute_engine
    *
    * For some reason this scope isn't listed as a constant under CloudKMSScopes.
    */
  val KmsScope = "https://www.googleapis.com/auth/cloudkms"

  /**
    * Scope to write metrics to Stackdriver Monitoring API.
    * Used by the monitoring action.
    *
    * For some reason we couldn't find this scope within Google libraries
    */
  val MonitoringWrite = "https://www.googleapis.com/auth/monitoring.write"

  /**
    * An image with the Google Cloud SDK installed.
    * http://gcr.io/google.com/cloudsdktool/cloud-sdk
    *
    * FYI additional older versions are available on DockerHub at:
    * https://hub.docker.com/r/google/cloud-sdk
    *
    * When updating this value, also consider updating the CromwellImagesSizeRoundedUpInGB below.
    */
  val CloudSdkImage = "gcr.io/google.com/cloudsdktool/cloud-sdk:276.0.0-slim"

  /*
   * At the moment, the cloud-sdk:slim (727MB on 2019-09-26) and possibly stedolan/jq (182MB) decompressed ~= 1 GB
   */
  val CromwellImagesSizeRoundedUpInGB = 1
}

case class ProjectLabels(labels: Map[String, String])
