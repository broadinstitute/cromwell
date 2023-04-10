package cromwell.backend.google.pipelines.batch.api

import com.google.cloud.batch.v1.AllocationPolicy.Accelerator
import com.google.cloud.batch.v1.{GetJobRequest, JobName}
import cromwell.backend.google.pipelines.batch.{BatchUtilityConversions, GcpBatchRequest, RunStatus}
import cromwell.core.WorkflowId
//import com.google.cloud.batch.v1.AllocationPolicy._
import com.google.cloud.batch.v1.AllocationPolicy.{InstancePolicy, InstancePolicyOrTemplate, LocationPolicy, NetworkInterface, NetworkPolicy, ProvisioningModel}
import com.google.cloud.batch.v1.{AllocationPolicy, ComputeResource, CreateJobRequest, Job, LogsPolicy, Runnable, ServiceAccount, TaskGroup, TaskSpec}
//import com.google.cloud.batch.v1.AllocationPolicy.{InstancePolicy, InstancePolicyOrTemplate, LocationPolicy, NetworkInterface, NetworkPolicy}
import com.google.cloud.batch.v1.AllocationPolicy.AttachedDisk
import com.google.cloud.batch.v1.LogsPolicy.Destination
import com.google.cloud.batch.v1.Runnable.Container
import com.google.cloud.batch.v1.Volume
import com.google.protobuf.Duration
import cromwell.backend.google.pipelines.batch.io.GcpBatchAttachedDisk
import org.slf4j.{Logger, LoggerFactory}

import scala.jdk.CollectionConverters._

class GcpBatchRequestFactoryImpl extends GcpBatchRequestFactory with BatchUtilityConversions {
  override def queryRequest(jobName: JobName): GetJobRequest = GetJobRequest.newBuilder.setName(jobName.toString).build

  val log: Logger = LoggerFactory.getLogger(RunStatus.toString)

  // VALUES HERE
  private val entryPoint = "/bin/sh"
  private val durationInSeconds: Long = 3600
  private val taskCount: Long = 1

  private def createRunnable(dockerImage: String, entryPoint: String, command: String): Runnable = {
    val runnable = Runnable.newBuilder.setContainer((Container.newBuilder.setImageUri(dockerImage).setEntrypoint(entryPoint).addCommands("-c").addCommands(command).build)).build

    runnable
  }

  private def createComputeResource(cpu: Long, memory: Long, bootDiskSizeMb: Long) = {
    ComputeResource
      .newBuilder
      .setCpuMilli(cpu)
      .setMemoryMib(memory)
      .setBootDiskMib(bootDiskSizeMb)
      .build
  }

  private def createInstancePolicy(machineType: String, cpuPlatform: String, spotModel: ProvisioningModel, accelerators: Option[Accelerator.Builder], attachedDisks: List[AttachedDisk]) = {

    //val attachedDiskIterable: Iterable[AttachedDisk] = attachedDisk.asJava
    //set GPU count to 0 if not included in workflow
    val gpuAccelerators = accelerators.getOrElse(Accelerator.newBuilder.setCount(0).setType(""))
    val instancePolicy = InstancePolicy
      .newBuilder
      .setMachineType(machineType)
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


  private def createNetworkInterface(vpcNetwork: String, vpcSubnetwork: String, noAddress: Boolean) = {
    NetworkInterface
      .newBuilder
      .setNoExternalIpAddress(noAddress)
      .setNetwork(vpcNetwork)
      .setSubnetwork(vpcSubnetwork)
      .build
  }

  private def createNetworkPolicy(networkInterface: NetworkInterface): NetworkPolicy = {
    NetworkPolicy
      .newBuilder
      .addNetworkInterfaces(0, networkInterface)
      .build
  }

  private def createTaskSpec(runnable: Runnable, computeResource: ComputeResource, retryCount: Int, durationInSeconds: Long, volumes: List[Volume]) = {
    TaskSpec
      .newBuilder
      .addRunnables(runnable)
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

  private def createAllocationPolicy(workflowId: WorkflowId, locationPolicy: LocationPolicy, instancePolicy: InstancePolicy, networkPolicy: NetworkPolicy, serviceAccount: ServiceAccount) = {
    AllocationPolicy
      .newBuilder
      .setLocation(locationPolicy)
      .setNetwork(networkPolicy)
      .putLabels("cromwell-workflow-id", workflowId.toString)
      .setServiceAccount(serviceAccount)
      .addInstances(InstancePolicyOrTemplate
        .newBuilder
        .setPolicy(instancePolicy)
        .build)
      .build
  }

  override def submitRequest(machineType: String, data: GcpBatchRequest): CreateJobRequest = {
    val batchAttributes = data.gcpBatchParameters.batchAttributes
    val runtimeAttributes = data.gcpBatchParameters.runtimeAttributes
    val createParameters = data.createParameters
    val retryCount = data.gcpBatchParameters.runtimeAttributes.preemptible
    val gcpBatchCommand: String = data.gcpBatchCommand
    val vpcNetwork: String = toVpcNetwork(batchAttributes)
    val vpcSubnetwork: String = toVpcSubnetwork(batchAttributes, runtimeAttributes)
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

    val noAddress = runtimeAttributes.noAddress

    println(f"command script container path ${data.createParameters.commandScriptContainerPath}")
    println(f"cloud workflow root ${data.createParameters.cloudWorkflowRoot}")
    println(f"all parameters ${data.createParameters.allParameters}")

    // parse preemption value and set value for Spot. Spot is replacement for preemptible
    val spotModel = toProvisioningModel(runtimeAttributes.preemptible)

    // Set GPU accelerators
    val accelerators = runtimeAttributes.gpuResource.map(toAccelerator)

    // Parse Service Account
    //val sa = batchAttributes.computeServiceAccount

    //val image = gcsTransferLibraryContainerPath
    //val gcsTransferLibraryContainerPath = createPipelineParameters.commandScriptContainerPath.sibling(GcsTransferLibraryName)
    val runnable = createRunnable(dockerImage = data.gcpBatchParameters.runtimeAttributes.dockerImage, entryPoint = entryPoint, command = gcpBatchCommand)

    val networkInterface = createNetworkInterface(vpcNetwork = vpcNetwork, vpcSubnetwork = vpcSubnetwork, noAddress = noAddress)
    val networkPolicy = createNetworkPolicy(networkInterface)
    val allDisks = toDisks(allDisksToBeMounted)
    val allVolumes = toVolumes(allDisksToBeMounted)
    val computeResource = createComputeResource(cpuCores, memory, gcpBootDiskSizeMb)
    val taskSpec = createTaskSpec(runnable, computeResource, retryCount, durationInSeconds, allVolumes)
    val taskGroup: TaskGroup = createTaskGroup(taskCount, taskSpec)
    val instancePolicy = createInstancePolicy(machineType = machineType, cpuPlatform = cpuPlatform, spotModel, accelerators, allDisks)
    val locationPolicy = LocationPolicy.newBuilder.addAllowedLocations(zones).build
    val allocationPolicy = createAllocationPolicy(data.workflowId, locationPolicy, instancePolicy.build, networkPolicy, gcpSa)
    val job = Job
      .newBuilder
      .addTaskGroups(taskGroup)
      .setAllocationPolicy(allocationPolicy)
      .putLabels("submitter", "cromwell") // label to signify job submitted by cromwell for larger tracking purposes within GCP batch
      .putLabels("cromwell-workflow-id", data.workflowId.toString) // label to make it easier to match Cromwell workflows with multiple GCP batch jobs
      .setLogsPolicy(LogsPolicy
        .newBuilder
        .setDestination(Destination
          .CLOUD_LOGGING)
        .build)

    CreateJobRequest
      .newBuilder
      .setParent(parent)
      .setJob(job)
      .setJobId(data.jobName)
      .build()
  }
}
