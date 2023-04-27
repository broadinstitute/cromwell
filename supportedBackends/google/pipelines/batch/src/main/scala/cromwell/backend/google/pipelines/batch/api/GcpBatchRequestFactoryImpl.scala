package cromwell.backend.google.pipelines.batch.api

import com.google.cloud.batch.v1.AllocationPolicy.Accelerator
import com.google.cloud.batch.v1.{DeleteJobRequest, GetJobRequest, JobName}
import cromwell.backend.google.pipelines.batch.models.GcpBatchConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.pipelines.batch.models.GcpBatchRequest
import cromwell.backend.google.pipelines.batch.runnable._
import cromwell.backend.google.pipelines.batch.util.BatchUtilityConversions
//import com.google.cloud.batch.v1.AllocationPolicy._
import com.google.cloud.batch.v1.AllocationPolicy.{AttachedDisk, InstancePolicy, InstancePolicyOrTemplate, LocationPolicy, NetworkInterface, NetworkPolicy, ProvisioningModel}
import com.google.cloud.batch.v1.LogsPolicy.Destination
import com.google.cloud.batch.v1.{AllocationPolicy, ComputeResource, CreateJobRequest, Job, LogsPolicy, Runnable, ServiceAccount, TaskGroup, TaskSpec, Volume}
import com.google.protobuf.Duration
import cromwell.backend.google.pipelines.batch.io.GcpBatchAttachedDisk
import cromwell.backend.google.pipelines.batch.models.VpcAndSubnetworkProjectLabelValues

import scala.jdk.CollectionConverters._

class GcpBatchRequestFactoryImpl()(implicit gcsTransferConfiguration: GcsTransferConfiguration) extends GcpBatchRequestFactory
  with BatchUtilityConversions
  with UserRunnable
  with ContainerSetup
  with Localization
  with Delocalization
  with MemoryRetryCheckRunnable
  with MonitoringRunnable
  with CheckpointingRunnable {

  override def queryRequest(jobName: JobName): GetJobRequest = GetJobRequest.newBuilder.setName(jobName.toString).build

  override def abortRequest(jobName: JobName): DeleteJobRequest = DeleteJobRequest.newBuilder.setName(jobName.toString).build()


  def createNetworkWithVPC(vpcAndSubnetworkProjectLabelValues: VpcAndSubnetworkProjectLabelValues, data: GcpBatchRequest): NetworkInterface.Builder = {

    val network = NetworkInterface
      .newBuilder
      .setNoExternalIpAddress(data.gcpBatchParameters.runtimeAttributes.noAddress)
      .setNetwork(vpcAndSubnetworkProjectLabelValues.networkName(data.gcpBatchParameters.projectId))

    vpcAndSubnetworkProjectLabelValues
      .subnetNameOption(data.gcpBatchParameters.projectId)
      .foreach(network.setSubnetwork)

    network

  }

  def createNetwork(data: GcpBatchRequest): NetworkInterface.Builder = {
    data.createParameters.vpcNetworkAndSubnetworkProjectLabels match {
      case Some(vpcAndSubnetworkProjectLabelValues) => createNetworkWithVPC(vpcAndSubnetworkProjectLabelValues, data)
      case _ => NetworkInterface.newBuilder().setNoExternalIpAddress(data.createParameters.runtimeAttributes.noAddress)
    }
  }

  private def createComputeResource(cpu: Long, memory: Long, bootDiskSizeMb: Long) = {
    ComputeResource
      .newBuilder
      .setCpuMilli(cpu)
      .setMemoryMib(memory)
      .setBootDiskMib(bootDiskSizeMb)
      .build
  }

  private def createInstancePolicy(cpuPlatform: String, spotModel: ProvisioningModel, accelerators: Option[Accelerator.Builder], attachedDisks: List[AttachedDisk]) = {

    //val attachedDiskIterable: Iterable[AttachedDisk] = attachedDisk.asJava

    //set GPU count to 0 if not included in workflow
    val gpuAccelerators = accelerators.getOrElse(Accelerator.newBuilder.setCount(0).setType(""))
    val instancePolicy = InstancePolicy
      .newBuilder
      .setProvisioningModel(spotModel)
      .addAllDisks(attachedDisks.asJava)
      .setMinCpuPlatform(cpuPlatform)
      .buildPartial()

    //add GPUs if GPU count is greater than 1
    if (gpuAccelerators.getCount >= 1) {
      val instancePolicyGpu = instancePolicy.toBuilder
      instancePolicyGpu.addAccelerators(gpuAccelerators).build
      instancePolicyGpu
    } else {
      instancePolicy.toBuilder
    }

  }


  private def createNetworkPolicy(networkInterface: NetworkInterface): NetworkPolicy = {
    NetworkPolicy
      .newBuilder
      .addNetworkInterfaces(0, networkInterface)
      .build
  }

  private def createTaskSpec(runnables: List[Runnable], computeResource: ComputeResource, retryCount: Int, durationInSeconds: Long, volumes: List[Volume]) = {
    TaskSpec
      .newBuilder
      .addAllRunnables(runnables.asJava)
      .setComputeResource(computeResource)
      .addAllVolumes(volumes.asJava)
      .setMaxRetryCount(retryCount)
      .setMaxRunDuration(Duration
        .newBuilder
        .setSeconds(durationInSeconds)
        .build)
  }

  private def createTaskGroup(taskCount: Long, task: TaskSpec.Builder): TaskGroup = {
    TaskGroup
      .newBuilder
      .setTaskCount(taskCount)
      .setTaskSpec(task)
      .build

  }

  private def createAllocationPolicy(data: GcpBatchRequest, locationPolicy: LocationPolicy, instancePolicy: InstancePolicy, networkPolicy: NetworkPolicy, serviceAccount: ServiceAccount) = {
    AllocationPolicy
      .newBuilder
      .setLocation(locationPolicy)
      .setNetwork(networkPolicy)
      .putLabels("cromwell-workflow-id", toLabel(data.workflowId.toString)) //label for workflow from WDL
      .putLabels("wdl-task-name", toLabel(data.gcpBatchParameters.jobDescriptor.taskCall.callable.name)) //label for task from WDL
      .putLabels("goog-batch-worker", "true")
      .setServiceAccount(serviceAccount)
      .addInstances(InstancePolicyOrTemplate
        .newBuilder
        .setPolicy(instancePolicy)
        .build)
      .build
  }

  override def submitRequest(data: GcpBatchRequest): CreateJobRequest = {

    val batchAttributes = data.gcpBatchParameters.batchAttributes
    val runtimeAttributes = data.gcpBatchParameters.runtimeAttributes
    val createParameters = data.createParameters
    val retryCount = data.gcpBatchParameters.runtimeAttributes.preemptible
    val allDisksToBeMounted: Seq[GcpBatchAttachedDisk] = createParameters.adjustedSizeDisks ++ createParameters.referenceDisksForLocalizationOpt.getOrElse(List.empty)
    val gcpBootDiskSizeMb = convertGbToMib(runtimeAttributes)

    // set parent for metadata storage of job information
    lazy val parent = s"projects/${data.gcpBatchParameters.projectId}/locations/${data.gcpBatchParameters.region}"
    val gcpSa = ServiceAccount.newBuilder.setEmail(batchAttributes.computeServiceAccount).build

    // make zones path
    val zones = toZonesPath(runtimeAttributes.zones)

    // convert to millicores for Batch
    val cpu = runtimeAttributes.cpu
    val cpuCores = toCpuCores(cpu.toString.toLong)

    val cpuPlatform = runtimeAttributes.cpuPlatform.getOrElse("")

    // convert memory to MiB for Batch
    val memory = toMemMib(runtimeAttributes.memory)

    // Determine max runtime for Batch
    val durationInSeconds: Long = data.gcpBatchParameters.batchAttributes.batchTimeout.toSeconds

    // Batch defaults to 1 task
    val taskCount: Long = 1

    println(f"command script container path ${data.createParameters.commandScriptContainerPath}")
    println(f"cloud workflow root ${data.createParameters.cloudWorkflowRoot}")
    println(f"all parameters:\n ${data.createParameters.allParameters.mkString("\n")}")

    // parse preemption value and set value for Spot. Spot is replacement for preemptible
    val spotModel = toProvisioningModel(runtimeAttributes.preemptible)

    // Set GPU accelerators
    val accelerators = runtimeAttributes.gpuResource.map(toAccelerator)
    
    //val image = gcsTransferLibraryContainerPath
    //val gcsTransferLibraryContainerPath = createPipelineParameters.commandScriptContainerPath.sibling(GcsTransferLibraryName)

    val networkInterface = createNetwork(data = data)
    val networkPolicy = createNetworkPolicy(networkInterface.build())
    val allDisks = toDisks(allDisksToBeMounted)
    val allVolumes = toVolumes(allDisksToBeMounted)

    val containerSetup: List[Runnable] = containerSetupRunnables(allVolumes)
    val localization: List[Runnable] = localizeRunnables(createParameters, allVolumes)
    val userRunnable: List[Runnable] = userRunnables(data.createParameters, allVolumes)
    val memoryRetryRunnable: List[Runnable] = checkForMemoryRetryRunnables(createParameters, allVolumes)
    val deLocalization: List[Runnable] = deLocalizeRunnables(createParameters, allVolumes)
    val monitoringSetup: List[Runnable] = monitoringSetupRunnables(createParameters, allVolumes)
    val monitoringShutdown: List[Runnable] = monitoringShutdownRunnables(createParameters)
    val checkpointingStart: List[Runnable] = checkpointingSetupRunnables(createParameters, allVolumes)
    val checkpointingShutdown: List[Runnable] = checkpointingShutdownRunnables(createParameters, allVolumes)
    val sshAccess: List[Runnable] = List.empty //sshAccessActions(createPipelineParameters, mounts)

    val sortedRunnables: List[Runnable] = RunnableUtils.sortRunnables(
      containerSetup = containerSetup,
        localization = localization,
        userRunnable = userRunnable,
        memoryRetryRunnable = memoryRetryRunnable,
        deLocalization = deLocalization,
        monitoringSetup = monitoringSetup,
        monitoringShutdown = monitoringShutdown,
        checkpointingStart = checkpointingStart,
        checkpointingShutdown = checkpointingShutdown,
        sshAccess = sshAccess,
        isBackground = _.getBackground,
      )

    val computeResource = createComputeResource(cpuCores, memory, gcpBootDiskSizeMb)
    val taskSpec = createTaskSpec(sortedRunnables, computeResource, retryCount, durationInSeconds, allVolumes)
    val taskGroup: TaskGroup = createTaskGroup(taskCount, taskSpec)
    val instancePolicy = createInstancePolicy(cpuPlatform, spotModel, accelerators, allDisks)
    val locationPolicy = LocationPolicy.newBuilder.addAllowedLocations(zones).build
    val allocationPolicy = createAllocationPolicy(data, locationPolicy, instancePolicy.build, networkPolicy, gcpSa)
    val job = Job
      .newBuilder
      .addTaskGroups(taskGroup)
      .setAllocationPolicy(allocationPolicy)
      .putLabels("submitter", "cromwell") // label to signify job submitted by cromwell for larger tracking purposes within GCP batch
      .putLabels("goog-batch-worker", "true")
      .putLabels("cromwell-workflow-id", toLabel(data.workflowId.toString)) // label to make it easier to match Cromwell workflows with multiple GCP batch jobs
      .putLabels("wdl-task-name", toLabel(data.gcpBatchParameters.jobDescriptor.taskCall.callable.name)) //label for task from WDL
      .setLogsPolicy(LogsPolicy
        .newBuilder
        .setDestination(Destination.CLOUD_LOGGING)
        .build)

    CreateJobRequest
      .newBuilder
      .setParent(parent)
      .setJob(job)
      .setJobId(data.jobName)
      .build()
  }
}
