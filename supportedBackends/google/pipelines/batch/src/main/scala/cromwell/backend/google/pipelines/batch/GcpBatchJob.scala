package cromwell.backend.google.pipelines.batch

import com.google.api.gax.rpc.{FixedHeaderProvider, HeaderProvider}
import com.google.cloud.batch.v1.AllocationPolicy.Accelerator
//import com.google.cloud.batch.v1.AllocationPolicy._
import com.google.cloud.batch.v1.AllocationPolicy.{InstancePolicy, InstancePolicyOrTemplate, LocationPolicy, NetworkInterface, NetworkPolicy, ProvisioningModel}
import com.google.cloud.batch.v1.{AllocationPolicy, BatchServiceClient, BatchServiceSettings, ComputeResource, CreateJobRequest, Job, LogsPolicy, Runnable, ServiceAccount, TaskGroup, TaskSpec}
//import com.google.cloud.batch.v1.AllocationPolicy.{InstancePolicy, InstancePolicyOrTemplate, LocationPolicy, NetworkInterface, NetworkPolicy}
import com.google.cloud.batch.v1.LogsPolicy.Destination
import com.google.cloud.batch.v1.Runnable.Container
import com.google.common.collect.ImmutableMap
import com.google.protobuf.Duration
import org.slf4j.{Logger, LoggerFactory}
//import scala.jdk.CollectionConverters._


final case class GcpBatchJob (
                             jobSubmission: GcpBatchRequest,
                             machineType: String
                            ) extends BatchUtilityConversions{

  val log: Logger = LoggerFactory.getLogger(RunStatus.toString)

  private val batchAttributes = jobSubmission.gcpBatchParameters.batchAttributes
  private val runtimeAttributes = jobSubmission.gcpBatchParameters.runtimeAttributes

  // VALUES HERE
  private val entryPoint = "/bin/sh"
  private val retryCount = jobSubmission.gcpBatchParameters.runtimeAttributes.preemptible
  private val durationInSeconds: Long = 3600
  private val taskCount: Long = 1
  private val gcpBatchCommand: String = jobSubmission.gcpBatchCommand
  private val vpcNetwork: String = toVpcNetwork(batchAttributes)
  private val vpcSubnetwork: String = toVpcSubnetwork(batchAttributes, runtimeAttributes)
  private val gcpBootDiskSizeMb = toBootDiskSizeMb(runtimeAttributes)

  // set user agent to cromwell so requests can be differentiated on batch
  private val user_agent_header = "user-agent"
  private val customUserAgentValue = "cromwell"
  private lazy val headerProvider: HeaderProvider = FixedHeaderProvider
    .create(ImmutableMap
      .of(user_agent_header, customUserAgentValue))

  private lazy val batchSettings = BatchServiceSettings.newBuilder.setHeaderProvider(headerProvider).build

  // TODO: Alex - Consider creating this client once, close it once this is not required
  lazy val batchServiceClient = BatchServiceClient.create(batchSettings)

  // set parent for metadata storage of job information
  lazy val parent = s"projects/${jobSubmission.gcpBatchParameters.projectId}/locations/${jobSubmission.gcpBatchParameters.region}"
  val gcpSa = ServiceAccount.newBuilder.setEmail(batchAttributes.computeServiceAccount).build

  // make zones path
  private val zones = toZonesPath(runtimeAttributes.zones)

  // convert to millicores for Batch
  private val cpu = runtimeAttributes.cpu
  private val cpuCores = toCpuCores(cpu.toString.toLong)

  private val cpuPlatform =  runtimeAttributes.cpuPlatform.getOrElse("")
  println(cpuPlatform)

  // convert memory to MiB for Batch
  private val memory = toMemMib(runtimeAttributes.memory)

  //private val bootDiskSize = runtimeAttributes.bootDiskSize
  private val noAddress = runtimeAttributes.noAddress

  println(f"command script container path ${jobSubmission.createParameters.commandScriptContainerPath}")
  println(f"cloud workflow root ${jobSubmission.createParameters.cloudWorkflowRoot}")
  println(f"all parameters ${jobSubmission.createParameters.allParameters}")

  // parse preemption value and set value for Spot. Spot is replacement for preemptible
  val spotModel = toProvisioningModel(runtimeAttributes.preemptible)

  // Set GPU accelerators
  private val accelerators = runtimeAttributes
    .gpuResource.map(toAccelerator)

  // Parse Service Account
  val sa = batchAttributes.computeServiceAccount

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

  private def createInstancePolicy(spotModel: ProvisioningModel, accelerators: Option[Accelerator.Builder]) = {

      //set GPU count to 0 if not included in workflow
      val gpuAccelerators = accelerators.getOrElse(Accelerator.newBuilder.setCount(0).setType(""))

      val instancePolicy = InstancePolicy
        .newBuilder
        .setMachineType(machineType)
        .setProvisioningModel(spotModel)
        //.setMinCpuPlatform(cpuPlatform)
        .buildPartial()

      //add GPUs if GPU count is greater than 1
      if(gpuAccelerators.getCount >= 1){
        val instancePolicyGpu = instancePolicy.toBuilder
        instancePolicyGpu.addAccelerators(gpuAccelerators).build
        instancePolicyGpu
      } else {
        instancePolicy.toBuilder
      }

  }


  private def createNetworkInterface(noAddress: Boolean) = {
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
      .build()
  }


  private def createTaskSpec(runnable: Runnable, computeResource: ComputeResource, retryCount: Int, durationInSeconds: Long) = {
    TaskSpec
      .newBuilder
      .addRunnables(runnable)
      .setComputeResource(computeResource)
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

  private def createAllocationPolicy(locationPolicy: LocationPolicy,  instancePolicy: InstancePolicy, networkPolicy: NetworkPolicy, serviceAccount: ServiceAccount) = {
    AllocationPolicy
      .newBuilder
      .setLocation(locationPolicy)
      .setNetwork(networkPolicy)
      .putLabels("cromwell-workflow-id", jobSubmission.workflowId.toString)
      .setServiceAccount(serviceAccount)
      .addInstances(InstancePolicyOrTemplate
        .newBuilder
        .setPolicy(instancePolicy)
        .build)
      .build
  }

  def submitJob(): Job = {
    //val image = gcsTransferLibraryContainerPath
    //val gcsTransferLibraryContainerPath = createPipelineParameters.commandScriptContainerPath.sibling(GcsTransferLibraryName)
    val runnable = createRunnable(dockerImage = jobSubmission.gcpBatchParameters.runtimeAttributes.dockerImage, entryPoint = entryPoint, command = gcpBatchCommand)

    val networkInterface = createNetworkInterface(noAddress)
    val networkPolicy = createNetworkPolicy(networkInterface)
    val computeResource = createComputeResource(cpuCores, memory, gcpBootDiskSizeMb)
    val taskSpec = createTaskSpec(runnable, computeResource, retryCount, durationInSeconds)
    val taskGroup: TaskGroup = createTaskGroup(taskCount, taskSpec)
    val instancePolicy = createInstancePolicy(spotModel, accelerators)
    val locationPolicy = LocationPolicy.newBuilder.addAllowedLocations(zones).build
    val allocationPolicy = createAllocationPolicy(locationPolicy, instancePolicy.build, networkPolicy, gcpSa)
    val job = Job
      .newBuilder
      .addTaskGroups(taskGroup)
      .setAllocationPolicy(allocationPolicy)
      .putLabels("submitter", "cromwell") // label to signify job submitted by cromwell for larger tracking purposes within GCP batch
      .putLabels("cromwell-workflow-id", jobSubmission.workflowId.toString) // label to make it easier to match Cromwell workflows with multiple GCP batch jobs
      .setLogsPolicy(LogsPolicy
        .newBuilder
        .setDestination(Destination
          .CLOUD_LOGGING)
        .build)

    val createJobRequest = CreateJobRequest
      .newBuilder
      .setParent(parent)
      .setJob(job)
      .setJobId(jobSubmission.jobName)
      .build()

    // TODO: Alex - Don't do blocking calls, consider using a Future instead
    val result = batchServiceClient
      .createJobCallable
      .call(createJobRequest)

    log.info("job submitted")
    batchServiceClient.close()
    log.info(result.getName)

    result.getUid
    result
  }
}