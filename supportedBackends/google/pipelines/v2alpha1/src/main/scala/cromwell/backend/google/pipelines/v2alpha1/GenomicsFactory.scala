package cromwell.backend.google.pipelines.v2alpha1

import java.net.URL

import com.google.api.client.http.{HttpRequest, HttpRequestInitializer}
import com.google.api.services.bigquery.BigqueryScopes
import com.google.api.services.compute.ComputeScopes
import com.google.api.services.genomics.v2alpha1.model._
import com.google.api.services.genomics.v2alpha1.{Genomics, GenomicsScopes}
import com.google.api.services.oauth2.Oauth2Scopes
import com.google.api.services.storage.StorageScopes
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.pipelines.common.action.ActionUtils
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.api.{PipelinesApiFactoryInterface, PipelinesApiRequestFactory}
import cromwell.backend.google.pipelines.common.{GoogleCloudScopes, VpcAndSubnetworkProjectLabelValues}
import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._
import cromwell.backend.google.pipelines.v2alpha1.api._
import cromwell.backend.standard.StandardAsyncJob
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.DockerConfiguration
import cromwell.core.logging.JobLogger
import mouse.all._
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize

import scala.jdk.CollectionConverters._

case class GenomicsFactory(applicationName: String, authMode: GoogleAuthMode, endpointUrl: URL)(implicit
  gcsTransferConfiguration: GcsTransferConfiguration
) extends PipelinesApiFactoryInterface
    with ContainerSetup
    with MonitoringAction
    with CheckpointingAction
    with Localization
    with UserAction
    with Delocalization
    with MemoryRetryCheckAction
    with SSHAccessAction {

  override def build(initializer: HttpRequestInitializer): PipelinesApiRequestFactory = new PipelinesApiRequestFactory {
    val genomics: Genomics = new Genomics.Builder(GoogleAuthMode.httpTransport, GoogleAuthMode.jsonFactory, initializer)
      .setApplicationName(applicationName)
      .setRootUrl(endpointUrl.toString)
      .build

    override def cancelRequest(job: StandardAsyncJob): HttpRequest =
      genomics.projects().operations().cancel(job.jobId, new CancelOperationRequest()).buildHttpRequest()

    override def getRequest(job: StandardAsyncJob): HttpRequest =
      genomics.projects().operations().get(job.jobId).buildHttpRequest()

    override def runRequest(createPipelineParameters: CreatePipelineParameters, jobLogger: JobLogger): HttpRequest = {
      def createNetworkWithVPC(vpcAndSubnetworkProjectLabelValues: VpcAndSubnetworkProjectLabelValues): Network = {
        val network = new Network()
          .setUsePrivateAddress(createPipelineParameters.runtimeAttributes.noAddress)
          .setName(vpcAndSubnetworkProjectLabelValues.networkName(createPipelineParameters.projectId))

        vpcAndSubnetworkProjectLabelValues
          .subnetNameOption(createPipelineParameters.projectId)
          .foreach(network.setSubnetwork)

        network
      }

      def createNetwork(): Network =
        createPipelineParameters.vpcNetworkAndSubnetworkProjectLabels match {
          case Some(vpcAndSubnetworkProjectLabelValues) => createNetworkWithVPC(vpcAndSubnetworkProjectLabelValues)
          case _ => new Network().setUsePrivateAddress(createPipelineParameters.runtimeAttributes.noAddress)
        }

      val allDisksToBeMounted = createPipelineParameters.disks ++
        createPipelineParameters.referenceDisksForLocalizationOpt.getOrElse(List.empty)

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
      val checkpointingStart: List[Action] = checkpointingSetupActions(createPipelineParameters, mounts)
      val checkpointingShutdown: List[Action] = checkpointingShutdownActions(createPipelineParameters)
      val sshAccess: List[Action] = sshAccessActions(createPipelineParameters, mounts)

      // adding memory as environment variables makes it easy for a user to retrieve the new value of memory
      // on the machine to utilize in their command blocks if needed
      val runtimeMemory = createPipelineParameters.runtimeAttributes.memory
      val environment =
        Map("MEM_UNIT" -> runtimeMemory.unit.toString, "MEM_SIZE" -> runtimeMemory.amount.toString).asJava

      val sortedActions =
        ActionUtils.sortActions[Action](
          containerSetup = containerSetup,
          localization = localization,
          userAction = userAction,
          memoryRetryAction = memoryRetryAction,
          deLocalization = deLocalization,
          monitoringSetup = monitoringSetup,
          monitoringShutdown = monitoringShutdown,
          checkpointingStart = checkpointingStart,
          checkpointingShutdown = checkpointingShutdown,
          sshAccess = sshAccess,
          isBackground = action =>
            Option(action.getFlags)
              .map(_.asScala)
              .toList
              .flatten
              .contains(ActionFlag.RunInBackground.toString)
        )

      val serviceAccount = new ServiceAccount()
        .setEmail(createPipelineParameters.computeServiceAccount)
        .setScopes(
          List(
            GenomicsScopes.GENOMICS,
            ComputeScopes.COMPUTE,
            StorageScopes.DEVSTORAGE_FULL_CONTROL,
            GoogleCloudScopes.KmsScope,
            // Profile and Email scopes are requirements for interacting with DRS Resolver
            Oauth2Scopes.USERINFO_EMAIL,
            Oauth2Scopes.USERINFO_PROFILE,
            // Monitoring scope as POC
            GoogleCloudScopes.MonitoringWrite,
            // Allow read/write with BigQuery
            BigqueryScopes.BIGQUERY
          ).asJava
        )

      val network: Network = createNetwork()

      val accelerators = createPipelineParameters.runtimeAttributes.gpuResource.map(toAccelerator).toList.asJava

      /*
       * Adjust using docker images used by Cromwell as well as the tool's docker image size if available
       */
      val adjustedBootDiskSize = {
        val fromRuntimeAttributes = createPipelineParameters.runtimeAttributes.bootDiskSize
        // Compute the decompressed size based on the information available
        val userCommandImageSizeInBytes = createPipelineParameters.jobDescriptor.dockerSize
          .map(_.toFullSize(DockerConfiguration.instance.sizeCompressionFactor))
          .getOrElse(0L)
        val userCommandImageSizeInGB =
          MemorySize(userCommandImageSizeInBytes.toDouble, MemoryUnit.Bytes).to(MemoryUnit.GB).amount
        val userCommandImageSizeRoundedUpInGB = userCommandImageSizeInGB.ceil.toInt

        val totalSize =
          fromRuntimeAttributes + userCommandImageSizeRoundedUpInGB + ActionUtils.cromwellImagesSizeRoundedUpInGB
        jobLogger.info(
          s"Adjusting boot disk size to $totalSize GB: $fromRuntimeAttributes GB (runtime attributes) + $userCommandImageSizeRoundedUpInGB GB (user command image) + ${ActionUtils.cromwellImagesSizeRoundedUpInGB} GB (Cromwell support images)"
        )
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

      if (createPipelineParameters.dockerImageCacheDiskOpt.isDefined) {
        jobLogger.info(
          "Docker image cache requested for the job, but the job is being executed by Google " +
            "Genomics API v2alpha1, while the feature is only supported by from Google Life Scielnces API starting from" +
            "the version v2beta"
        )
      }

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
        .setTimeout(createPipelineParameters.pipelineTimeout.toSeconds.toString + "s")

      val pipelineRequest = new RunPipelineRequest()
        .setPipeline(pipeline)
        .setLabels(createPipelineParameters.googleLabels.map(label => label.key -> label.value).toMap.asJava)

      genomics.pipelines().run(pipelineRequest).buildHttpRequest()
    }
  }

  override def usesEncryptedDocker: Boolean = true
}
