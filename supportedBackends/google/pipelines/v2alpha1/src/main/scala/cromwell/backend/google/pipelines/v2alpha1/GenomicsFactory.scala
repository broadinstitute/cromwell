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
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.{LocalizationConfiguration, VirtualPrivateCloudConfiguration}
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.api.{PipelinesApiFactoryInterface, PipelinesApiRequestFactory}
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

case class GenomicsFactory(applicationName: String, authMode: GoogleAuthMode, endpointUrl: URL)(implicit localizationConfiguration: LocalizationConfiguration) extends PipelinesApiFactoryInterface
  with ContainerSetup
  with MonitoringAction
  with Localization
  with UserAction
  with Delocalization {

  override def build(initializer: HttpRequestInitializer): PipelinesApiRequestFactory = new PipelinesApiRequestFactory {
    implicit lazy val googleProjectMetadataLabelDecoder: Decoder[ProjectLabels] = deriveDecoder

    val ResourceManagerAuthScopes = List(GenomicsScopes.CLOUD_PLATFORM)
    val VirtualPrivateCloudNetworkPath = "projects/%s/global/networks/%s/"

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
              .setUsePrivateAddress(createPipelineParameters.runtimeAttributes.noAddress)
              .setName(VirtualPrivateCloudNetworkPath.format(createPipelineParameters.projectId, networkLabel._2))

            subnetworkLabelOption foreach { case(_, subnet) => network.setSubnetwork(subnet) }
            network
          case None =>
            // Falling back to running the job on default network since the project does not provide the custom network
            // specifying keys in its metadata
            new Network().setUsePrivateAddress(createPipelineParameters.runtimeAttributes.noAddress)
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
          case None => new Network().setUsePrivateAddress(createPipelineParameters.runtimeAttributes.noAddress)
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
      val deLocalization: List[Action] = deLocalizeActions(createPipelineParameters, mounts)
      val monitoring: List[Action] = monitoringActions(createPipelineParameters, mounts)
      val allActions = containerSetup ++ localization ++ userAction ++ deLocalization ++ monitoring

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
            Oauth2Scopes.USERINFO_PROFILE,
            // Monitoring scope as POC
            GenomicsFactory.MonitoringWrite,
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
        /*
         * At the moment, google/cloud-sdk:slim (664MB on 12/3/18) and possibly stedolan/jq (182MB) decompressed ~= 1 GB
         */
        val cromwellImagesSize = 1
        val fromRuntimeAttributes = createPipelineParameters.runtimeAttributes.bootDiskSize
        // Compute the decompressed size based on the information available
        val actualSizeInBytes = createPipelineParameters.jobDescriptor.dockerSize.map(_.toFullSize(DockerConfiguration.instance.sizeCompressionFactor)).getOrElse(0L)
        val actualSizeInGB = MemorySize(actualSizeInBytes.toDouble, MemoryUnit.Bytes).to(MemoryUnit.GB).amount
        val actualSizeRoundedUpInGB = actualSizeInGB.ceil.toInt

        // If the size of the image is larger than what is in the runtime attributes, adjust it to uncompressed size
        if (actualSizeRoundedUpInGB > fromRuntimeAttributes) {
          jobLogger.info(s"Adjusting boot disk size to $actualSizeRoundedUpInGB GB to account for docker image size")
        }

        Math.max(fromRuntimeAttributes, actualSizeRoundedUpInGB) + cromwellImagesSize
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

      val pipelineRequest = new RunPipelineRequest()
        .setPipeline(pipeline)
        .setLabels(createPipelineParameters.googleLabels.map(label => label.key -> label.value).toMap.asJava)

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

  /**
    * Scope to write metrics to Stackdriver Monitoring API.
    * Used by the monitoring action.
    *
    * For some reason we couldn't find this scope within Google libraries
    */
  val MonitoringWrite = "https://www.googleapis.com/auth/monitoring.write"
}

case class ProjectLabels(labels: Map[String, String])
