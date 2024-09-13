package cromwell.backend.google.batch.api

import com.google.cloud.batch.v1.AllocationPolicy._
import com.google.cloud.batch.v1.LogsPolicy.Destination
import com.google.cloud.batch.v1.{
  AllocationPolicy,
  ComputeResource,
  CreateJobRequest,
  DeleteJobRequest,
  GetJobRequest,
  Job,
  JobName,
  LogsPolicy,
  Runnable,
  ServiceAccount,
  TaskGroup,
  TaskSpec,
  Volume
}
import com.google.protobuf.Duration
import cromwell.backend.google.batch.io.GcpBatchAttachedDisk
import cromwell.backend.google.batch.models.GcpBatchConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.batch.models.{GcpBatchLogsPolicy, GcpBatchRequest, VpcAndSubnetworkProjectLabelValues}
import cromwell.backend.google.batch.runnable._
import cromwell.backend.google.batch.util.{BatchUtilityConversions, GcpBatchMachineConstraints}
import cromwell.core.labels.{Label, Labels}
import cromwell.core.logging.JobLogger

import scala.jdk.CollectionConverters._

class GcpBatchRequestFactoryImpl()(implicit gcsTransferConfiguration: GcsTransferConfiguration)
    extends GcpBatchRequestFactory
    with BatchUtilityConversions
    with UserRunnable
    with ContainerSetup
    with Localization
    with Delocalization
    with MemoryRetryCheckRunnable
    with MonitoringRunnable
    with CheckpointingRunnable {

  override def queryRequest(jobName: JobName): GetJobRequest = GetJobRequest.newBuilder.setName(jobName.toString).build

  override def abortRequest(jobName: JobName): DeleteJobRequest =
    DeleteJobRequest.newBuilder.setName(jobName.toString).build()

  def createNetworkWithVPC(vpcAndSubnetworkProjectLabelValues: VpcAndSubnetworkProjectLabelValues,
                           data: GcpBatchRequest
  ): NetworkInterface.Builder = {

    val network = NetworkInterface.newBuilder
      .setNoExternalIpAddress(data.gcpBatchParameters.runtimeAttributes.noAddress)
      .setNetwork(vpcAndSubnetworkProjectLabelValues.networkName(data.gcpBatchParameters.projectId))

    vpcAndSubnetworkProjectLabelValues
      .subnetNameOption(data.gcpBatchParameters.projectId)
      .foreach(network.setSubnetwork)

    network

  }

  def createNetwork(data: GcpBatchRequest): NetworkInterface.Builder =
    data.createParameters.vpcNetworkAndSubnetworkProjectLabels match {
      case Some(vpcAndSubnetworkProjectLabelValues) => createNetworkWithVPC(vpcAndSubnetworkProjectLabelValues, data)
      case _ => NetworkInterface.newBuilder().setNoExternalIpAddress(data.createParameters.runtimeAttributes.noAddress)
    }

  private def createComputeResource(cpu: Long, memory: Long, bootDiskSizeMb: Long) =
    ComputeResource.newBuilder
      .setCpuMilli(cpu)
      .setMemoryMib(memory)
      .setBootDiskMib(bootDiskSizeMb)
      .build

  private def createInstancePolicy(cpuPlatform: String,
                                   spotModel: ProvisioningModel,
                                   accelerators: Option[Accelerator.Builder],
                                   attachedDisks: List[AttachedDisk],
                                   machineType: String
  ): InstancePolicy.Builder = {

    // set GPU count to 0 if not included in workflow
    val gpuAccelerators = accelerators.getOrElse(Accelerator.newBuilder.setCount(0).setType("")) // TODO: Driver version

    val instancePolicy = InstancePolicy.newBuilder
      .setProvisioningModel(spotModel)
      .setMachineType(machineType)
      .addAllDisks(attachedDisks.asJava)
      .setMinCpuPlatform(cpuPlatform)
      .buildPartial()

    // add GPUs if GPU count is greater than 1
    if (gpuAccelerators.getCount >= 1) {
      val instancePolicyGpu = instancePolicy.toBuilder
      instancePolicyGpu.addAccelerators(gpuAccelerators).build
      instancePolicyGpu
    } else {
      instancePolicy.toBuilder
    }

  }

  private def createNetworkPolicy(networkInterface: NetworkInterface): NetworkPolicy =
    NetworkPolicy.newBuilder
      .addNetworkInterfaces(0, networkInterface)
      .build

  private def createTaskSpec(runnables: List[Runnable],
                             computeResource: ComputeResource,
                             retryCount: Int,
                             durationInSeconds: Long,
                             volumes: List[Volume]
  ) =
    TaskSpec.newBuilder
      .addAllRunnables(runnables.asJava)
      .setComputeResource(computeResource)
      .addAllVolumes(volumes.asJava)
      .setMaxRetryCount(retryCount)
      .setMaxRunDuration(
        Duration.newBuilder
          .setSeconds(durationInSeconds)
          .build
      )

  private def createTaskGroup(taskCount: Long, task: TaskSpec.Builder): TaskGroup =
    TaskGroup.newBuilder
      .setTaskCount(taskCount)
      .setTaskSpec(task)
      .build

  private def createAllocationPolicy(data: GcpBatchRequest,
                                     locationPolicy: LocationPolicy,
                                     instancePolicy: InstancePolicy,
                                     networkPolicy: NetworkPolicy,
                                     serviceAccount: ServiceAccount,
                                     accelerators: Option[Accelerator.Builder]
  ) = {

    val allocationPolicy = AllocationPolicy.newBuilder
      .setLocation(locationPolicy)
      .setNetwork(networkPolicy)
      .putLabels("cromwell-workflow-id", toLabel(data.workflowId.toString)) // label for workflow from WDL
      .putLabels("goog-batch-worker", "true")
      .putAllLabels(data.createParameters.googleLabels.map(label => label.key -> label.value).toMap.asJava)
      .setServiceAccount(serviceAccount)
      .buildPartial()

    val gpuAccelerators = accelerators.getOrElse(Accelerator.newBuilder.setCount(0).setType(""))

    // add GPUs if GPU count is greater than or equal to 1
    if (gpuAccelerators.getCount >= 1) {
      allocationPolicy.toBuilder.addInstances(
        InstancePolicyOrTemplate.newBuilder.setPolicy(instancePolicy).setInstallGpuDrivers(true).build
      )
    } else {
      allocationPolicy.toBuilder.addInstances(InstancePolicyOrTemplate.newBuilder.setPolicy(instancePolicy).build)
    }
  }

  override def submitRequest(data: GcpBatchRequest, jobLogger: JobLogger): CreateJobRequest = {

    val runtimeAttributes = data.gcpBatchParameters.runtimeAttributes
    val createParameters = data.createParameters
    val retryCount = data.gcpBatchParameters.runtimeAttributes.preemptible
    val allDisksToBeMounted: Seq[GcpBatchAttachedDisk] =
      createParameters.disks ++ createParameters.referenceDisksForLocalizationOpt.getOrElse(List.empty)
    val gcpBootDiskSizeMb = convertGbToMib(runtimeAttributes)

    // set parent for metadata storage of job information
    lazy val parent = s"projects/${createParameters.projectId}/locations/${data.gcpBatchParameters.region}"
    val gcpSa = ServiceAccount.newBuilder.setEmail(createParameters.computeServiceAccount).build

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

    // parse preemption value and set value for Spot. Spot is replacement for preemptible
    val spotModel = toProvisioningModel(runtimeAttributes.preemptible)

    // Set GPU accelerators
    val accelerators = runtimeAttributes.gpuResource.map(toAccelerator)

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
    val sshAccess: List[Runnable] = List.empty // sshAccessActions(createPipelineParameters, mounts)

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
      isBackground = _.getBackground
    )

    val computeResource = createComputeResource(cpuCores, memory, gcpBootDiskSizeMb)
    val taskSpec = createTaskSpec(sortedRunnables, computeResource, retryCount, durationInSeconds, allVolumes)
    val taskGroup: TaskGroup = createTaskGroup(taskCount, taskSpec)
    val machineType = GcpBatchMachineConstraints.machineType(runtimeAttributes.memory,
                                                             runtimeAttributes.cpu,
                                                             cpuPlatformOption = runtimeAttributes.cpuPlatform,
                                                             jobLogger = jobLogger
    )
    val instancePolicy =
      createInstancePolicy(cpuPlatform = cpuPlatform, spotModel, accelerators, allDisks, machineType = machineType)
    val locationPolicy = LocationPolicy.newBuilder.addAllowedLocations(zones).build
    val allocationPolicy =
      createAllocationPolicy(data, locationPolicy, instancePolicy.build, networkPolicy, gcpSa, accelerators)
    val logsPolicy = data.gcpBatchParameters.batchAttributes.logsPolicy match {
      case GcpBatchLogsPolicy.CloudLogging =>
        LogsPolicy.newBuilder.setDestination(Destination.CLOUD_LOGGING).build
      case GcpBatchLogsPolicy.Path =>
        ???
    }

    val googleLabels = data.createParameters.googleLabels.map(l => Label(l.key, l.value))

    val jobDescriptor = data.createParameters.jobDescriptor
    val backendJobDescriptorKey = jobDescriptor.key

    val workflow = jobDescriptor.workflowDescriptor
    val call = jobDescriptor.taskCall
    val subWorkflow = workflow.callable
    val subWorkflowLabels =
      if (!subWorkflow.equals(workflow.rootWorkflow))
        Labels("cromwell-sub-workflow-name" -> subWorkflow.name,
               "cromwell-sub-workflow-id" -> s"cromwell-sub-${jobDescriptor.workflowDescriptor.id.toString}"
        )
      else
        Labels.empty

    val alias = call.localName
    val aliasLabels =
      if (!alias.equals(call.callable.name))
        Labels("wdl-call-alias" -> alias)
      else
        Labels.empty

    val shardLabels = Labels(backendJobDescriptorKey.index.map(l => Label("wdl-shard-index", l.toString)).toVector)

    val allLabels = Labels(
      "cromwell-workflow-id" -> s"cromwell-${workflow.rootWorkflowId}",
      "wdl-task-name" -> call.callable.name,
      "wdl-attempt" -> backendJobDescriptorKey.attempt.toString,
      "goog-batch-worker" -> "true",
      "submitter" -> "cromwell"
    ) ++ shardLabels ++ subWorkflowLabels ++ aliasLabels ++ Labels(googleLabels.toVector)

    val job = Job.newBuilder
      .addTaskGroups(taskGroup)
      .setAllocationPolicy(allocationPolicy.build())
      .putAllLabels(allLabels.asJavaMap)
      .setLogsPolicy(logsPolicy)

    CreateJobRequest.newBuilder
      .setParent(parent)
      .setJob(job)
      .setJobId(data.jobName)
      .build()
  }
}
